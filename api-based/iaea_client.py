"""Client for the IAEA Live Chart of Nuclides REST API (nds.iaea.org).

Provides the two datasets the API-based tools need:

- ``ground_states()``: one row per nuclide -- half-life, decay modes and
  branching ratios.
- ``build_gamma_library()``: gamma lines (energy, intensity, parent state,
  half-life) for every parent nuclide in a configurable half-life window,
  assembled into the record format the spectrum analyzer matches against.

Every response is cached on disk under ``cache/``, so the full library build
costs one pass of ~1000 small requests once and is free afterwards. The
assembled library is also written to ``data/api_gamma_library.csv``.
"""

import io
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
CACHE = HERE / 'cache'
DATA = HERE / 'data'

BASE_URL = 'https://nds.iaea.org/relnsd/v1/data?'
USER_AGENT = 'Mozilla/5.0 (X11; Linux x86_64) nuclear-data-client'
REQUEST_PAUSE_S = 0.12          # polite pacing between uncached requests
RETRIES = 3

# Library selection: parent states and gamma lines worth matching against.
MIN_PARENT_HALF_LIFE_S = 60.0
MAX_PARENT_HALF_LIFE_S = 1E17     # keeps primordial naturals such as 40K
MIN_INTENSITY_PERCENT = 0.01
ENERGY_RANGE_KEV = (40.0, 2000.0)
MAX_MASS_NUMBER = 236


def _fetch(query):
    """GET one API query, serving from the disk cache when possible."""
    CACHE.mkdir(exist_ok=True)
    cache_file = CACHE / (query.replace('&', '_').replace('=', '-') + '.csv')
    if cache_file.exists():
        return cache_file.read_text(encoding='utf-8')

    request = urllib.request.Request(BASE_URL + query,
                                     headers={'User-Agent': USER_AGENT})
    last_error = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                text = response.read().decode('utf-8', 'replace')
            cache_file.write_text(text, encoding='utf-8')
            time.sleep(REQUEST_PAUSE_S)
            return text
        except (urllib.error.URLError, OSError) as error:
            last_error = error
            time.sleep(2.0 * (attempt + 1))
    raise ConnectionError(f'IAEA API request failed after {RETRIES} tries: '
                          f'{query} ({last_error})')


def ground_states():
    """All nuclides with half-life, decay modes and branching ratios."""
    frame = pd.read_csv(io.StringIO(_fetch('fields=ground_states&nuclides=all')))
    frame['half_life_sec'] = pd.to_numeric(frame['half_life_sec'],
                                           errors='coerce')
    return frame


def decay_gammas(nuclide):
    """Gamma lines from the decay of every state of one nuclide, e.g. '99tc'."""
    text = _fetch(f'fields=decay_rads&nuclides={nuclide}&rad_types=g')
    if not text.strip() or text.splitlines()[0].startswith('<'):
        return pd.DataFrame()
    try:
        return pd.read_csv(io.StringIO(text))
    except pd.errors.ParserError:
        return pd.DataFrame()


def _half_life_units(seconds):
    """Split a half-life in seconds into the analyzer's (value, unit) form."""
    for limit, unit, divisor in [(60.0, 's', 1.0),
                                 (3600.0, 'm', 60.0),
                                 (86400.0, 'h', 3600.0),
                                 (31557600.0, 'd', 86400.0)]:
        if seconds < limit:
            return seconds / divisor, unit
    return seconds / 31557600.0, 'y'


def candidate_parents(states=None):
    """Nuclides whose ground state can survive long enough to be measured."""
    states = ground_states() if states is None else states
    hl = states['half_life_sec']
    mass = states['z'] + states['n']
    selected = states[(hl >= MIN_PARENT_HALF_LIFE_S)
                      & (hl <= MAX_PARENT_HALF_LIFE_S)
                      & (mass <= MAX_MASS_NUMBER)
                      & (states['z'] >= 3)]
    return [f'{int(row.z + row.n)}{row.symbol.lower()}'
            for row in selected.itertuples()]


def build_gamma_library(progress=None):
    """Assemble the gamma library from per-nuclide decay radiation data.

    Returns a DataFrame with columns ``energy`` (keV), ``intensity`` (%),
    ``isotope`` (e.g. ``99Tcm``), ``half_life`` and ``unit`` -- the same record
    shape as the curated text library. Cached to ``data/api_gamma_library.csv``.
    """
    DATA.mkdir(exist_ok=True)
    library_file = DATA / 'api_gamma_library.csv'
    if library_file.exists():
        return pd.read_csv(library_file)

    low_kev, high_kev = ENERGY_RANGE_KEV
    records = []
    parents = candidate_parents()
    for index, nuclide in enumerate(parents):
        if progress and index % 50 == 0:
            progress(index, len(parents))
        gammas = decay_gammas(nuclide)
        if gammas.empty or 'energy' not in gammas.columns:
            continue
        for column in ('energy', 'intensity', 'p_energy', 'half_life_sec'):
            gammas[column] = pd.to_numeric(gammas[column], errors='coerce')
        gammas = gammas.dropna(subset=['energy', 'intensity', 'half_life_sec'])
        gammas = gammas[(gammas.energy >= low_kev) & (gammas.energy <= high_kev)
                        & (gammas.intensity >= MIN_INTENSITY_PERCENT)
                        & (gammas.half_life_sec >= MIN_PARENT_HALF_LIFE_S)
                        & (gammas.half_life_sec <= MAX_PARENT_HALF_LIFE_S)]
        for row in gammas.itertuples():
            value, unit = _half_life_units(row.half_life_sec)
            name = f'{int(row.p_z + row.p_n)}{row.p_symbol}'
            if row.p_energy > 0:
                name += 'm'         # line comes from a metastable parent state
            records.append((round(float(row.energy), 4),
                            round(float(row.intensity), 6),
                            name, round(value, 4), unit))

    library = pd.DataFrame(records, columns=['energy', 'intensity', 'isotope',
                                             'half_life', 'unit'])
    library = library.drop_duplicates().sort_values('energy')
    library.to_csv(library_file, index=False)
    return library


def library_records(library=None):
    """The library as the list-of-lists records the analyzer matches against."""
    library = build_gamma_library() if library is None else library
    return [[float(row.energy), float(row.intensity), str(row.isotope),
             float(row.half_life), str(row.unit)]
            for row in library.itertuples()]
