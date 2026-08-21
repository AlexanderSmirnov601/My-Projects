# My-Projects

Three tools written for nuclear-physics work: an ETL job that turns the NNDC
NuDat export into a normalized MySQL database, a desktop viewer for the chart of
nuclides built on top of it, and a gamma-spectrum analyzer.

| Project | What it does |
| --- | --- |
| [Online mySQL application](<Online mySQL application>) | Parses the NNDC/NuDat CSV export, derives the Z/N change of every decay mode, and loads it into the `IsotopeDB` schema. |
| [Radioisotope map & Activity calculation](<Radioisotope map & Activity calculation>) | Tkinter app: half-life heatmap of the chart of nuclides, decay-chain diagrams and Bateman activity curves. |
| [Spectrum Analyzer](<Spectrum Analyzer>) | Reads Ortec `.Spe` spectra, finds and fits peaks, identifies isotopes against a gamma library, and computes activity and production yield. |

## Setup

```bash
python -m venv .venv
.venv/Scripts/activate        # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

On Debian/Ubuntu the Tkinter viewer also needs the `python3-tk` system package.

## Configuration

No credentials or machine-specific paths are stored in the source. Each project
reads a config file next to its script — `settings.ini` for the Spectrum
Analyzer, `config.ini` for the two database projects — and each supports a
gitignored `*.local.ini` override for machine-specific values.

Database passwords are never written to a file. They are read from
`ISOTOPEDB_PASSWORD` and `ISOTOPEDB_GUEST_PASSWORD`, and prompted for when unset.
