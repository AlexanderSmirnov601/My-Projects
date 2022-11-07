import mysql.connector
import pandas as pd
import numpy as np
import getpass
lala = []


def zn_change(line):
    z = np.nan
    n = np.nan
    if 'F' in line or 'STABLE' in line:
        return z, n
    elif "IT" in line:
        z = 0
        n = 0
    elif 'EC' in line:
        try:
            multiplier = int(line[0])
            z = -2
            n = 2
        except:
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
            try:
                multiplier = int(letter)
            except:
                if letter == 'A':
                    z += multiplier*(-2)
                    n += multiplier*(-2)
                    multiplier = 1
                if letter == 'P':
                    z += multiplier*(-1)
                    multiplier = 1
                if letter == 'N':
                    n += multiplier*(-1)
                    multiplier = 1
                if letter == 'B':
                    z += multiplier*(1)
                    n += multiplier*(-1)
                    multiplier = 1
                if letter == 'E':
                    z += multiplier*(-1)
                    n += multiplier*(1)
                    multiplier = 1

    return z, n

DB_password = getpass.getpass('Enter a database password: ')

nudat = mysql.connector.connect(
    host='34.67.94.150',
    user='root',
    password = DB_password,
    )

cursor = nudat.cursor()
cursor.execute('''
                DROP DATABASE IF EXISTS IsotopeDB;
                ''')
cursor.execute('''
                CREATE DATABASE IsotopeDB;
                ''')
cursor.execute('''
                USE IsotopeDB;
                ''')
cursor.execute('''
               DROP USER IF EXISTS 'guest'@'%';
               ''')
cursor.execute('''
               CREATE USER 'guest'@'%'
               IDENTIFIED BY 'password';
               ''')
cursor.execute('''
               GRANT SELECT ON * TO 'guest'@'%'; 
               ''')
cursor.execute('''
                FLUSH PRIVILEGES;
                ''')

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
nudat.commit()


wanted_columns = ['halflife', 'decayModes', 'z', 'name', 'n', 'levelEnergy(MeV)']
og_data = pd.read_csv('./nndc_nudat_data_export.csv')
col = og_data.columns
flag = False
for column in og_data.columns:
    for check in wanted_columns:
        if column == check:
            flag = True
    if flag:
        flag = False
        continue
    else:
        og_data.drop(column, axis=1, inplace=True)


og_data['dm'] = og_data['decayModes'].str[0:2]

irrelevant_dms = ['Mg', 'Ne']
time = [['as', 1E-18],
        ['fs', 1E-15],
        ['ps', 1E-12],
        ['ns', 1E-9],
        ['us', 1E-6],
        ['ms', 1E-3],
        ['s', 1],
        ['m', 60],
        ['h', 3600],
        ['d', 24*3600],
        ['y', 365.25*24*3600]]

for dm in irrelevant_dms:
    og_data = og_data.drop(og_data[og_data.dm == dm].index)


data = og_data.loc[:, 'decayModes']
data = data.dropna()
data = data.sort_values()


unique_DMs = ['STABLE']
for line in data.unique():
    if pd.isna(line):
        continue
    unique_DMs.append(line.split(' ')[0])
unique_DMs = np.unique(unique_DMs)
unique_DMs = np.sort(unique_DMs)
for dm in unique_DMs:
    dm = str(dm)
    if np.isnan(zn_change(dm)[0]):
        cursor.execute('''INSERT INTO DecayTypes(name)
                       VALUES (%s);''', (dm,))
    else:
        cursor.execute('''INSERT INTO DecayTypes(name, z_change, n_change)
                       VALUES (%s, %s, %s);''', (dm, *zn_change(dm)))


cursor.execute('''SELECT id, name, z_change, n_change FROM DecayTypes;''')
decayModes = cursor.fetchall()


for index, row in og_data.iterrows():
    z = row['z']
    n = row['n']
    name = row['name']
    hl_calc = []
    dec_prob = None
    dm_id = None
    halflife = None
    child_z = None
    child_n = None
    e_level_mev = None


    if not pd.isna(row['decayModes']) and '?' in row['decayModes'].split(' '):
        og_data = og_data.drop([index])
        continue
    elif not pd.isna(row['decayModes']) and pd.isna(row['halflife']):
        og_data = og_data.drop([index])
        continue
    if pd.isna(row['decayModes']):
        og_data['decayModes'][index] = 'STABLE'

    for dm in decayModes:
        if og_data['decayModes'][index].split(' ')[0] == dm[1]:
            dm_id = dm[0]
            if not pd.isna(dm[2]):
                child_z = z + int(dm[2])
                child_n = n + int(dm[3])
            break


    try:
        dec_prob = float(row['decayModes'].split(' ')[-1])
    except:
        dec_prob = None

    if pd.isna(row['halflife']):
        og_data = og_data.drop([index])
        continue
    elif row['halflife'] == 'STABLE':
        halflife = 'STABLE'
    else:
        hl_calc = row['halflife'].split(' ')
        # print(hl_calc)
        try:
            hl_calc[0] = float(hl_calc[0])
        except:
            hl_calc.pop(0)
            hl_calc[0] = float(hl_calc[0])
        if hl_calc[1] in ('mev', 'kev', 'ev'):
            og_data = og_data.drop([index])
            continue
        for t in time:
            if hl_calc[1] == t[0]:
                halflife = hl_calc[0]*t[1]
                break
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
                       (hl_sec, probability, decaytype_id, isotope_z, isotope_n, e_level_mev)
                       VALUES (NULL, NULL, %s, %s, %s, %s);''',
                       (dm_id, z, n, e_level_mev)
                       )
    elif dec_prob is None:
        cursor.execute('''INSERT INTO DecayModes
                       (hl_sec, probability, decaytype_id, isotope_z, isotope_n, child_z, child_n, e_level_mev)
                       VALUES (%s, NULL, %s, %s, %s, %s, %s, %s);''',
                       (halflife, dm_id, z, n, child_z, child_n, e_level_mev)
                       )
    else:
        cursor.execute('''INSERT INTO DecayModes
                       (hl_sec, probability, decaytype_id, isotope_z, isotope_n, child_z, child_n, e_level_mev)
                       VALUES (%s, %s, %s, %s, %s, %s, %s, %s);''',
                       (halflife, dec_prob, dm_id, z, n, child_z, child_n, e_level_mev)
                       )


    nudat.commit()
cursor.close()
nudat.close()