# My-Projects

Three tools from nuclear-physics work: a builder that turns the NNDC NuDat
export into a normalized MySQL database, a desktop viewer for the chart of
nuclides built on that database, and a gamma-spectrum analyzer for activation
measurements.

| Project | What it does |
| --- | --- |
| [isotope-database-builder](isotope-database-builder) | Parses the NNDC/NuDat CSV export, derives the Z/N change of every decay mode, and loads the result into the `IsotopeDB` MySQL schema. |
| [nuclide-chart-viewer](nuclide-chart-viewer) | Tkinter app over `IsotopeDB`: half-life heatmap of the chart of nuclides, decay-chain diagrams, and activity-vs-time curves for a chosen isotope. |
| [gamma-spectrum-analyzer](gamma-spectrum-analyzer) | Reads Ortec `.Spe` spectra, finds and fits peaks, identifies isotopes against a gamma library, and computes activity and production yield. |
| [api-based](api-based) | Both tools above rebuilt on the public IAEA Live Chart of Nuclides API — no local database or curated library needed. |

Each of the two graphical projects ships a `sample_outputs.pdf` collecting its
figures with short explanations.

## Setup

```bash
python -m venv .venv
.venv/Scripts/activate        # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

On Debian/Ubuntu the Tkinter viewer also needs the `python3-tk` system package.

## Configuration

No credentials or machine-specific paths live in the source. Each project reads
a config file next to its script — `settings.ini` for the spectrum analyzer,
`config.ini` for the two database projects — and each supports a gitignored
`*.local.ini` override for machine-specific values.

Database passwords are never written to a file: they are read from
`ISOTOPEDB_PASSWORD` / `ISOTOPEDB_GUEST_PASSWORD`, with an interactive prompt
as fallback.
