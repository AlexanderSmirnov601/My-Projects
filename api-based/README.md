# api-based

The gamma-spectrum analyzer and the nuclide chart viewer rebuilt on the
public [IAEA Live Chart of Nuclides API](https://nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html)
(`nds.iaea.org`) instead of local data files and a private MySQL database.

Both tools import their pipeline and figure code from the original projects —
only the data layer differs — so results are directly comparable between the
file-based and API-based versions.

| Script | Replaces | Data source |
| --- | --- | --- |
| `analyze_spectrum_api.py` | `gamma_library.txt` | per-nuclide decay radiation (`fields=decay_rads`) |
| `nuclide_chart_api.py` | `IsotopeDB` (MySQL) | bulk ground states (`fields=ground_states&nuclides=all`) |

## Usage

```bash
python analyze_spectrum_api.py
```

```bash
python nuclide_chart_api.py
```

No database or credentials are needed. Spectrum settings (file, peaks,
sample parameters) still come from `../gamma-spectrum-analyzer/settings.ini`.

## Data handling

- Every API response is cached under `cache/` (gitignored), so repeated runs
  are offline and instant.
- The assembled gamma library is committed as `data/api_gamma_library.csv`.
  Delete it to force a rebuild from the live API (~1000 small requests,
  a few minutes). Selection windows (parent half-life, line intensity,
  energy range) are constants at the top of `iaea_client.py`.

## Differences from the file-based versions

The curated `gamma_library.txt` carries ~2k hand-picked lines; the API library
carries every catalogued line in the selection window (~50k), so ambiguous
peaks can resolve differently. Two systematic effects are worth knowing:

- **Equilibrium half-lives.** The curated library was built by hand with
  transient equilibrium deliberately encoded in the half-life field: a
  generator-fed daughter is listed under the half-life of the nuclide that
  rate-limits its decay in a real sample, not its own. Its 140.5 keV entry
  carries the 66 h half-life of the feeding 99Mo rather than the 6 h of 99mTc
  itself, which is exactly how that activity behaves in a sample days after
  irradiation — so the line still ranks correctly when the raw state
  half-life would have written it off. The API reports each state's own
  half-life, and the ranking then discards short-lived daughters that are in
  fact still present through their parent. Reproducing this properly in the
  API version would mean modelling parent feeding in the ranking; it would
  arrive at the same answer the curation encodes in a text file.
- **Natural background.** The ranking scores lines by production-and-decay
  plausibility, which correctly favours activation products — and therefore
  ranks primordial background lines (40K at 1461 keV) below them, where the
  small curated library had no competitors in the window.

Nuclide data otherwise reflects current IAEA evaluations rather than the
2013-era NuDat export; several half-lives that were placeholders there are
measured values here.
