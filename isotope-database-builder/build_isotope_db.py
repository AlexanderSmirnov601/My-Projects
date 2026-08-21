"""Build the IsotopeDB MySQL database from an NNDC/NuDat CSV export.

Creates the DecayTypes / Isotope / DecayModes schema, derives the Z and N change
of every decay mode so child nuclides can be resolved, converts half-lives to
seconds and loads the lot.

Connection details live in ``config.ini`` next to this file; passwords are read
from the environment (``ISOTOPEDB_PASSWORD``, ``ISOTOPEDB_GUEST_PASSWORD``) and
prompted for when unset.
"""

import configparser
import getpass
import os
import re
from pathlib import Path

import mysql.connector
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
IDENTIFIER = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')

# Marks a "half-life" quoted as a level width in eV; those rows are skipped.
LEVEL_WIDTH = object()

# Listed in NuDat export order so the loaded frame matches the source layout.
WANTED_COLUMNS = ['z', 'n', 'name', 'levelEnergy(MeV)', 'halflife',
                  'decayModes']
IRRELEVANT_DECAY_MODES = ['Mg', 'Ne']
TIME_UNITS = [['as', 1E-18],
              ['fs', 1E-15],
              ['ps', 1E-12],
              ['ns', 1E-9],
              ['us', 1E-6],
              ['ms', 1E-3],
              ['s', 1],
              ['m', 60],
              ['h', 3600],
              ['d', 24 * 3600],
              ['y', 365.25 * 24 * 3600]]


def load_config():
    config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    path = HERE / 'config.ini'
    if not config.read([path, path.with_suffix('.local.ini')],
                       encoding='utf-8'):
        raise FileNotFoundError(f'no config file found at {path}')
    return config


def secret(env_var, prompt):
    """Password from the environment, falling back to an interactive prompt."""
    return os.environ.get(env_var) or getpass.getpass(prompt)


def zn_change(line):
    """Change in proton and neutron number produced by a decay mode string.

    ``"B-"`` -> (1, -1), ``"2A"`` -> (-4, -4), ``"24Mg"`` -> (-12, -12), and
    ``(nan, nan)`` for modes such as spontaneous fission where no single child
    nuclide can be derived.
    """
    z = np.nan
    n = np.nan
    if 'F' in line or 'STABLE' in line:
        return z, n
    elif "IT" in line:
        z = 0
        n = 0
    elif 'EC' in line:
        # A leading digit marks double electron capture, e.g. "2EC".
        if line[0].isdigit():
            z = -2
            n = 2
        else:
            z = -1
            n = 1
    elif 'Ne' in line:
        z = - 10
        n = - (int(line[0:2]) + z)
    elif 'Mg' in line:
        z = - 12
        n = - (int(line[0:2]) + z)
    elif 'C' in line:
        z = - 6
        n = - (int(line[0:2]) + z)
    elif 'O' in line:
        z = - 8
        n = - (int(line[0:2]) + z)
    elif 'Si' in line:
        z = - 14
        n = - (int(line[0:2]) + z)
    else:
        z = 0
        n = 0
        multiplier = 1
        for letter in line:
            if letter.isdigit():
                multiplier = int(letter)
                continue
            if letter == 'A':
                z += multiplier * (-2)
                n += multiplier * (-2)
                multiplier = 1
            if letter == 'P':
                z += multiplier * (-1)
                multiplier = 1
            if letter == 'N':
                n += multiplier * (-1)
                multiplier = 1
            if letter == 'B':
                z += multiplier * (1)
                n += multiplier * (-1)
                multiplier = 1
            if letter == 'E':
                z += multiplier * (-1)
                n += multiplier * (1)
                multiplier = 1

    return z, n


def halflife_seconds(halflife):
    """Convert an NNDC half-life string to seconds.

    Returns ``'STABLE'`` for stable nuclides, a float for a recognised time
    unit, ``LEVEL_WIDTH`` for entries quoted as a level width in eV (those rows
    are dropped) and ``None`` for an unrecognised unit (stored as NULL).
    """
    if halflife == 'STABLE':
        return 'STABLE'

    parts = halflife.split(' ')
    try:
        value = float(parts[0])
    except ValueError:
        parts.pop(0)                    # leading comparator such as ">" or "~"
        value = float(parts[0])

    if parts[1] in ('mev', 'kev', 'ev'):
        return LEVEL_WIDTH

    for unit, factor in TIME_UNITS:
        if parts[1] == unit:
            return value * factor
    return None


def create_schema(cursor, database, guest_user, guest_password):
    if not IDENTIFIER.match(database):
        raise ValueError(f'unsafe database name in config.ini: {database!r}')

    cursor.execute(f'DROP DATABASE IF EXISTS {database};')
    cursor.execute(f'CREATE DATABASE {database};')
    cursor.execute(f'USE {database};')
    cursor.execute('DROP USER IF EXISTS %s@\'%\';', (guest_user,))
    cursor.execute('CREATE USER %s@\'%\' IDENTIFIED BY %s;',
                   (guest_user, guest_password))
    cursor.execute('GRANT SELECT ON * TO %s@\'%\';', (guest_user,))
    cursor.execute('FLUSH PRIVILEGES;')

    cursor.execute('''
        CREATE TABLE DecayTypes(
            id INTEGER auto_increment NOT NULL UNIQUE,
            name TEXT,
            z_change INTEGER,
            n_change INTEGER,
            PRIMARY KEY(id)
            )
        ''')
    cursor.execute('''
        CREATE TABLE Isotope(
            z INTEGER,
            n INTEGER,
            name TEXT,
            CONSTRAINT id PRIMARY KEY (z, n)
            );
        ''')
    cursor.execute('''
        CREATE TABLE DecayModes(
            id INTEGER auto_increment NOT NULL UNIQUE,
            hl_sec DOUBLE,
            probability FLOAT,
            decaytype_id INTEGER,
            isotope_z INTEGER NOT NULL,
            isotope_n INTEGER NOT NULL,
            child_z INTEGER,
            child_n INTEGER,
            e_level_mev FLOAT,
            FOREIGN KEY(decaytype_id) REFERENCES DecayTypes(id),
            CONSTRAINT isotope_id FOREIGN KEY (isotope_z, isotope_n)
                            REFERENCES Isotope(z, n),
            PRIMARY KEY(id)
            );
        ''')


def load_export(path):
    """Read the NuDat export, keeping only the columns we actually use."""
    og_data = pd.read_csv(path)[WANTED_COLUMNS]
    og_data['dm'] = og_data['decayModes'].str[0:2]
    for dm in IRRELEVANT_DECAY_MODES:
        og_data = og_data.drop(og_data[og_data.dm == dm].index)
    return og_data


def insert_decay_types(cursor, og_data):
    """One DecayTypes row per distinct decay mode, with its Z/N change."""
    data = og_data.loc[:, 'decayModes'].dropna().sort_values()

    unique_dms = ['STABLE']
    for line in data.unique():
        if pd.isna(line):
            continue
        unique_dms.append(line.split(' ')[0])
    unique_dms = np.sort(np.unique(unique_dms))

    for dm in unique_dms:
        dm = str(dm)
        if np.isnan(zn_change(dm)[0]):
            cursor.execute('''INSERT INTO DecayTypes(name)
                           VALUES (%s);''', (dm,))
        else:
            cursor.execute('''INSERT INTO DecayTypes(name, z_change, n_change)
                           VALUES (%s, %s, %s);''', (dm, *zn_change(dm)))

    cursor.execute('''SELECT id, name, z_change, n_change FROM DecayTypes;''')
    return cursor.fetchall()


def insert_isotopes(cursor, connection, og_data, decay_types):
    for index, row in og_data.iterrows():
        z = row['z']
        n = row['n']
        name = row['name']
        dm_id = None
        child_z = None
        child_n = None

        if not pd.isna(row['decayModes']) and '?' in row['decayModes'].split(' '):
            og_data = og_data.drop([index])
            continue
        elif not pd.isna(row['decayModes']) and pd.isna(row['halflife']):
            og_data = og_data.drop([index])
            continue

        decay_modes = row['decayModes']
        if pd.isna(decay_modes):
            # A missing decayModes entry marks a stable nuclide.
            decay_modes = 'STABLE'
            og_data.loc[index, 'decayModes'] = decay_modes

        for dm in decay_types:
            if decay_modes.split(' ')[0] == dm[1]:
                dm_id = dm[0]
                if not pd.isna(dm[2]):
                    child_z = z + int(dm[2])
                    child_n = n + int(dm[3])
                break

        try:
            dec_prob = float(row['decayModes'].split(' ')[-1])
        except (ValueError, AttributeError):
            dec_prob = None

        if pd.isna(row['halflife']):
            og_data = og_data.drop([index])
            continue

        halflife = halflife_seconds(row['halflife'])
        if halflife is LEVEL_WIDTH:
            og_data = og_data.drop([index])
            continue

        e_level_mev = row['levelEnergy(MeV)']
        if np.isnan(e_level_mev):
            e_level_mev = 0

        cursor.execute('''INSERT IGNORE INTO Isotope
                       (z, n, name)
                       VALUES (%s, %s, %s);''',
                       (z, n, name)
                       )
        if halflife == 'STABLE':
            cursor.execute('''INSERT INTO DecayModes
                           (hl_sec, probability, decaytype_id, isotope_z,
                            isotope_n, e_level_mev)
                           VALUES (NULL, NULL, %s, %s, %s, %s);''',
                           (dm_id, z, n, e_level_mev)
                           )
        else:
            cursor.execute('''INSERT INTO DecayModes
                           (hl_sec, probability, decaytype_id, isotope_z,
                            isotope_n, child_z, child_n, e_level_mev)
                           VALUES (%s, %s, %s, %s, %s, %s, %s, %s);''',
                           (halflife, dec_prob, dm_id, z, n, child_z, child_n,
                            e_level_mev)
                           )
        connection.commit()


def main():
    config = load_config()
    database = config['database']['database']

    export_path = Path(config['paths']['nudat_export'])
    if not export_path.is_absolute():
        export_path = HERE / export_path

    connection = mysql.connector.connect(
        host=config['database']['host'],
        user=config['database']['user'],
        password=secret('ISOTOPEDB_PASSWORD', 'Enter a database password: '),
        )
    try:
        cursor = connection.cursor()
        try:
            create_schema(
                cursor, database, config['guest']['user'],
                secret('ISOTOPEDB_GUEST_PASSWORD',
                       'Enter a password for the read-only guest account: '))
            connection.commit()

            og_data = load_export(export_path)
            decay_types = insert_decay_types(cursor, og_data)
            insert_isotopes(cursor, connection, og_data, decay_types)
            connection.commit()
            print(f'{database} rebuilt from {export_path.name}')
        finally:
            cursor.close()
    finally:
        connection.close()


if __name__ == '__main__':
    main()
