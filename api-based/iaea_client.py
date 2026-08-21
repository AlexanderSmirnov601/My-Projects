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

# Post-pull curation: equilibrium half-life propagation. A daughter fed by a
# longer-lived parent decays with the parent's half-life in a real sample
# (transient equilibrium), so each line is assigned the half-life of its
# rate-limiting upstream feeder. Feeding is followed through branches of at
# least EQUILIBRIUM_MIN_BRANCH_PERCENT, and a feeder only propagates its
# half-life if it is at most EQUILIBRIUM_MAX_PARENT_S -- longer-lived
# ancestors (primordial chains) do not drive activation-sample equilibria.
EQUILIBRIUM_MIN_BRANCH_PERCENT = 1.0
EQUILIBRIUM_MAX_PARENT_S = 3.156E7          # one year

# Post-pull curation: ambient background lines. These are always present in a
# counting room regardless of the sample and do not decay on measurement
# timescales, so they are pinned to a fixed effective half-life that keeps
# them competitive in the ranking -- the curated library achieved the same by
# carrying them with no competitors in the match window.
BACKGROUND_EFFECTIVE_HALF_LIFE = (30.0, 'd')
KNOWN_BACKGROUND_LINES = [                  # (energy keV, isotope)
    (1460.82, '40K'),
    (661.657, '137Cs'),
    (609.32, '214Bi'), (1120.29, '214Bi'), (1764.49, '214Bi'),
    (351.93, '214Pb'), (295.22, '214Pb'),
    (583.19, '208Tl'), (911.2, '228Ac'), (968.97, '228Ac'),
    (186.21, '226Ra'), (238.63, '212Pb'),
]

# Change in (Z, N) produced by each IAEA decay-mode label; None marks modes
# with no single child nuclide (fission).
MODE_CHANGES = {
    'B-': (1, -1), 'B+': (-1, 1), 'EC': (-1, 1), 'EC+B+': (-1, 1),
    'A': (-2, -2), 'IT': (0, 0),
    'B-N': (1, -2), 'B-2N': (1, -3), 'B-A': (-1, -3),
    'ECP': (-2, 1), 'B+P': (-2, 1),
    'N': (0, -1), '2N': (0, -2), 'P': (-1, 0), '2P': (-2, 0),
    '2B-': (2, -2), '2EC': (-2, 2), '2B+': (-2, 2),
    'SF': None, 'ECSF': None,
}


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
            name = f'{int(row.p_z + row.p_n)}{row.p_symbol}'
            if row.p_energy > 0:
                name += 'm'         # line comes from a metastable parent state
            records.append((round(float(row.energy), 4),
                            round(float(row.intensity), 6),
                            name, int(row.p_z), int(row.p_n),
                            float(row.half_life_sec)))

    library = pd.DataFrame(records, columns=['energy', 'intensity', 'isotope',
                                             'z', 'n', 'state_hl_sec'])
    library = library.drop_duplicates().sort_values('energy')
    library = curate_library(library)
    library.to_csv(library_file, index=False)
    return library


def feeder_half_lives(states=None):
    """Rate-limiting upstream feeder half-life per nuclide, in seconds.

    Builds the nuclide-level feeding graph from the ground-states decay modes
    and propagates each feeder's decay timescale down its chains to a fixed
    point, so a chain such as 140Ba -> 140La carries 140Ba's 12.75 d to the
    140La lines just as 99Mo carries its 66 h to the 99mTc line.
    """
    states = ground_states() if states is None else states
    states = states.rename(columns={'decay_1_%': 'decay_1_pct',
                                    'decay_2_%': 'decay_2_pct',
                                    'decay_3_%': 'decay_3_pct'})

    own = {}
    children = {}
    for row in states.itertuples():
        if str(row.half_life) == 'STABLE' or pd.isna(row.half_life_sec):
            continue
        z, n = int(row.z), int(row.n)
        own[(z, n)] = float(row.half_life_sec)
        for mode, percent in [(row.decay_1, row.decay_1_pct),
                              (row.decay_2, row.decay_2_pct),
                              (row.decay_3, row.decay_3_pct)]:
            change = MODE_CHANGES.get(str(mode))
            if (change is None or change == (0, 0) or pd.isna(percent)
                    or float(percent) < EQUILIBRIUM_MIN_BRANCH_PERCENT):
                continue
            children.setdefault((z, n), []).append((z + change[0],
                                                    n + change[1]))

    # up[X]: timescale on which X's activity disappears from a sample --
    # its own half-life or that of its slowest capped feeder.
    up = dict(own)
    for _ in range(60):
        changed = False
        for parent, kids in children.items():
            drive = up[parent]
            if drive > EQUILIBRIUM_MAX_PARENT_S:
                continue
            for kid in kids:
                if kid in up and drive > up[kid]:
                    up[kid] = drive
                    changed = True
        if not changed:
            break

    feeders = {}
    for parent, kids in children.items():
        drive = up[parent]
        if drive > EQUILIBRIUM_MAX_PARENT_S:
            continue
        for kid in kids:
            feeders[kid] = max(feeders.get(kid, 0.0), drive)
    return feeders


def curate_library(library, states=None):
    """The post-pull stage: equilibrium half-lives and background pinning.

    Adds ``eff_hl_sec`` (the half-life the ranking should use), the analyzer's
    ``half_life``/``unit`` split of it, and a ``background`` flag. The raw
    ``state_hl_sec`` column is kept so the manipulation stays inspectable.
    """
    feeders = feeder_half_lives(states)
    effective = [max(row.state_hl_sec, feeders.get((row.z, row.n), 0.0))
                 for row in library.itertuples()]
    library = library.assign(eff_hl_sec=effective, background=False)

    value, unit = BACKGROUND_EFFECTIVE_HALF_LIFE
    seconds = value * {'s': 1, 'm': 60, 'h': 3600, 'd': 86400,
                       'y': 31557600}[unit]
    for energy, isotope in KNOWN_BACKGROUND_LINES:
        mask = ((library.isotope == isotope)
                & ((library.energy - energy).abs() < 0.5))
        library.loc[mask, 'eff_hl_sec'] = seconds
        library.loc[mask, 'background'] = True

    split = [_half_life_units(s) for s in library.eff_hl_sec]
    library = library.assign(
        half_life=[round(v, 4) for v, _ in split],
        unit=[u for _, u in split])
    return library


def library_records(library=None):
    """The library as the list-of-lists records the analyzer matches against."""
    library = build_gamma_library() if library is None else library
    return [[float(row.energy), float(row.intensity), str(row.isotope),
             float(row.half_life), str(row.unit)]
            for row in library.itertuples()]
