# isotope-database-builder

Turns the NNDC/NuDat CSV export into the normalized `IsotopeDB` MySQL database
that [nuclide-chart-viewer](../nuclide-chart-viewer) reads.

## Schema

| Table | Contents |
| --- | --- |
| `DecayTypes` | One row per distinct decay mode, with the change in proton and neutron number it produces. |
| `Isotope` | Nuclides, keyed on `(z, n)`. |
| `DecayModes` | Half-life in seconds, branching ratio, excitation level, and the resulting child nuclide, per isotope. |

`zn_change()` parses NNDC decay-mode strings (`B-`, `2A`, `EC`, `24Mg`) into
the Z and N deltas the decay produces, so child nuclides can be resolved and
decay chains walked. Modes with no single child nuclide, such as spontaneous
fission, store `NULL`.

## Running

```bash
python build_isotope_db.py
```

This **drops and recreates** the database — a full rebuild, not an incremental
load. It also creates the read-only `guest` account the viewer uses.

Set `ISOTOPEDB_PASSWORD` (admin) and `ISOTOPEDB_GUEST_PASSWORD` (the account to
create), or enter both at the prompt. Host, user, database name and the CSV
path are set in `config.ini`; the gitignored `config.local.ini` overrides it.

## Data

`nndc_nudat_data_export.csv` is an export from the NNDC NuDat database. Six of
its 38 columns are used: `z`, `n`, `name`, `levelEnergy(MeV)`, `halflife` and
`decayModes`.
