# Online mySQL application

Turns the NNDC/NuDat CSV export into a normalized MySQL database that the
[Radioisotope map viewer](<../Radioisotope map & Activity calculation>) reads.

## Schema

| Table | Contents |
| --- | --- |
| `DecayTypes` | One row per distinct decay mode, with the change in proton and neutron number it produces. |
| `Isotope` | Nuclides, keyed on `(z, n)`. |
| `DecayModes` | Half-life in seconds, branching ratio, excitation level and the resulting child nuclide, per isotope. |

The interesting part is `zn_change()`, which parses NNDC decay-mode strings
(`B-`, `2A`, `EC`, `24Mg`) into the Z and N deltas the decay produces, so child
nuclides can be resolved and chains walked. Modes with no single child nuclide,
such as spontaneous fission, store `NULL`.

## Running

```bash
python "CSV converter.py"
```

This **drops and recreates** the database, so it is a rebuild rather than an
incremental load. It also creates the read-only `guest` account the viewer uses.

Set `ISOTOPEDB_PASSWORD` (the admin account) and `ISOTOPEDB_GUEST_PASSWORD` (the
account to be created), or enter both at the prompt. Host, user, database name
and the CSV path are set in `config.ini`; `config.local.ini` overrides it and is
gitignored.

## Data

`nndc_nudat_data_export.csv` is an export from the NNDC NuDat database. Only six
of its 38 columns are used: `z`, `n`, `name`, `levelEnergy(MeV)`, `halflife` and
`decayModes`.
