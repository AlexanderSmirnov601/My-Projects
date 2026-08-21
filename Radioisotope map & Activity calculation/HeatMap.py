"""Tkinter viewer for the IsotopeDB nuclide chart.

Three tabs: a half-life heatmap of the whole chart of nuclides with a zoomed
inset on the isotope you enter, the decay chain leading from it to a stable
nuclide, and the Bateman activity of every member of that chain over time.

Connection details live in ``config.ini`` next to this file; the read-only
account's password comes from ``ISOTOPEDB_GUEST_PASSWORD`` or an interactive
prompt.
"""

import configparser
import getpass
import os
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

import matplotlib.pyplot as plt
import mysql.connector
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from PIL import Image, ImageTk
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent

# Chart of nuclides extent: 120 proton rows by 180 neutron columns.
CHART_PROTONS = 120
CHART_NEUTRONS = 180
STABLE_PLACEHOLDER = 1E+30      # stand-in half-life so stable nuclides colour in
MAX_CHAIN_LENGTH = 200          # guards against cycles through isomeric states


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #
def load_config():
    config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    path = HERE / 'config.ini'
    if not config.read([path, path.with_suffix('.local.ini')],
                       encoding='utf-8'):
        raise FileNotFoundError(f'no config file found at {path}')
    return config


def use_style(name, fallback='seaborn-v0_8'):
    """Apply a matplotlib style, tolerating the 3.6 seaborn-style rename."""
    for candidate in (name, fallback):
        try:
            plt.style.use(candidate)
            return
        except OSError:
            continue
    print(f'warning: neither {name!r} nor {fallback!r} is an available style')


# --------------------------------------------------------------------------- #
# database
# --------------------------------------------------------------------------- #
def isotope_chain_data(z, n, data, cursor, visited=None):
    """Walk the decay chain from (z, n) down to a stable nuclide.

    ``visited`` stops the walk from looping forever when a chain cycles back on
    itself through an isomeric transition; without it a cycle in the data means
    unbounded recursion.
    """
    if visited is None:
        visited = set()
    if (z, n) in visited or len(visited) >= MAX_CHAIN_LENGTH:
        return data
    visited.add((z, n))

    cursor.execute(
                    '''SELECT m.hl_sec, m.isotope_z,
                    m.isotope_n, i.name, t.name,
                    m.child_z, m.child_n, m.probability,
                    m.e_level_mev
                    FROM DecayModes m join Isotope i  join DecayTypes t
                    where m.isotope_z = i.z
                    and m.isotope_n = i.n
                    and t.id = m.decaytype_id
                    and m.isotope_z=%s and m.isotope_n=%s
                    order by m.e_level_mev asc, m.probability desc;''', (z, n)
                    )
    arr = cursor.fetchall()
    if not arr:
        return data

    if arr[0][5]:
        data.append(arr[0])
        isotope_chain_data(arr[0][5], arr[0][6], data, cursor, visited)
    elif arr[0][4] == 'STABLE' and len(arr) >= 2 and arr[1][4] == 'IT':
        data.append(arr[1])
        data.append(arr[0])
        return data
    else:
        data.append(arr[0])
        return data

    if len(arr) >= 2 and arr[1][7] > 1 and arr[1][-1] == 0:
        data.append(arr[1])
        isotope_chain_data(arr[1][5], arr[1][6], data, cursor, visited)
    return data


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def heatmap_img(data, xy, results_dir):
    """Half-life heatmap of the whole chart, with an inset around ``xy``."""
    arr = np.full((CHART_PROTONS, CHART_NEUTRONS), np.nan)
    annotation = [[None for _ in range(CHART_NEUTRONS)]
                  for _ in range(CHART_PROTONS)]

    for row in data:
        half_life, z, n, isotope, decay_type = row[0], row[1], row[2], row[3], row[4]
        arr[z][n] = STABLE_PLACEHOLDER if half_life is None else half_life
        annotation[z][n] = f'{isotope}\n{decay_type}'
        if half_life is not None:
            annotation[z][n] += '\n' + '{:.2e}'.format(half_life) + ' s'

    log_norm = LogNorm(vmin=1E-3, vmax=1E+15)

    use_style('seaborn')
    fig, _ = plt.subplots(figsize=[21, 15])

    ax = sns.heatmap(arr, cmap='YlGnBu_r',
                     xticklabels=10, yticklabels=10,
                     norm=log_norm, square=True,
                     cbar_kws={"shrink": 0.76, "pad": 0.01}
                     )
    ax.invert_yaxis()
    ax.set_facecolor("white")
    ax.set_xlabel('Number of neutrons, N =', fontsize=20,
                  fontname='Cambria'
                  )
    ax.set_ylabel('Number of protons, Z =', fontsize=20,
                  fontname='Cambria'
                  )
    plt.tick_params(axis='x', which='both', top=False)
    plt.tick_params(axis='y', which='both', right=False)
    plt.xticks(fontsize=12, fontname='Cambria')
    plt.yticks(fontsize=12, fontname='Cambria')

    c_bar = ax.collections[0].colorbar
    c_bar.set_label('Half-life time, seconds', fontsize=20,
                    fontname='Cambria'
                    )
    c_bar.ax.tick_params(labelsize=12)
    for labels in c_bar.ax.yaxis.get_ticklabels():
        labels.set_family("Cambria")

    plt.xlim([1, CHART_NEUTRONS - 2])
    plt.ylim([0, CHART_PROTONS - 2])

    axins = ax.inset_axes([0.595, 0.058, 0.45, 0.53])
    for z, line in enumerate(annotation):
        for n, _ in enumerate(line):
            if z < xy[1] - 2 or z > xy[1] + 2 or n < xy[0] - 2 or n > xy[0] + 2:
                annotation[z][n] = np.nan
                arr[z][n] = np.nan

    ax1 = sns.heatmap(arr, ax=axins, cmap='YlGnBu_r',
                      xticklabels=1, yticklabels=1,
                      linecolor='black', linewidths=2,
                      cbar=False, fmt='',
                      annot=annotation,
                      annot_kws=dict(clip_on=True,
                                     weight='bold',
                                     fontname='Cambria'
                                     ),
                      norm=log_norm, square=True
                      )
    ax1.invert_yaxis()
    ax1.set_facecolor("white")

    axins.set_xlim(max(xy[0] - 2, 0), xy[0] + 3)
    axins.set_ylim(max(xy[1] - 2, 0), xy[1] + 3)

    axins.set_xticklabels(axins.get_xmajorticklabels(), fontsize=14,
                          weight='bold', fontname='Cambria'
                          )
    axins.set_yticklabels(axins.get_ymajorticklabels(), fontsize=14,
                          weight='bold', fontname='Cambria'
                          )

    _, pp1, pp2 = mark_inset(ax, axins,
                             loc1=1, loc2=1, edgecolor='firebrick',
                             linewidth=3, alpha=1
                             )
    pp1.loc1 = 2
    pp2.loc1 = 2
    pp1.loc2 = pp2.loc2 = 1 if xy[0] < 120 else 3

    path = results_dir / 'HeatMap.png'
    plt.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def chain_img(data, line_counter, x, y, ax, depth=0):
    """Draw the decay chain as a left-to-right ladder of boxed nuclides."""
    if depth >= MAX_CHAIN_LENGTH or line_counter >= len(data):
        return

    ax.add_patch(Rectangle((x, y), 1, 1, fill=True,
                           ec='black', fc='skyblue', linewidth=2
                           )
                 )
    ax.text(x + .5, y + .5, data[line_counter][3],
            fontsize=15, fontname='Cambria',
            color="black", ha="center",
            va="center", weight='bold'
            )

    # enumerate, not data.index(line): index() returns the FIRST row equal to
    # this one, which is the wrong branch whenever a chain repeats a nuclide.
    for offset, line in enumerate(data[line_counter + 1:], line_counter + 1):
        if (
            line[3] == data[line_counter][3]
            and data[offset - 1][4] == 'STABLE'
        ):
            ax.annotate("", xy=(x + .5, y - 1), xytext=(x + .5, y - .25),
                        arrowprops=dict(width=3, fc='black'),
                        fontname='Cambria'
                        )
            ax.text(x + .25, y - .625,
                    str(line[7]) + '%',
                    fontsize=15, color="black", ha="center",
                    va="center", rotation=-90, fontname='Cambria'
                    )
            ax.text(x + .7, y - .625,
                    "{:.1e}".format(line[0]) + ' s',
                    fontsize=15, color="black", ha="center",
                    va="center", rotation=-90, fontname='Cambria'
                    )
            ax.text(x + .95, y - .625, str(line[4]) + ' decay',
                    color="black", ha="center", va="center", style='italic',
                    rotation=-90, fontsize=15, fontname='Cambria'
                    )
            plt.plot(x, y - 2.25)
            chain_img(data, offset + 1, x, y - 2.25, ax, depth + 1)
            break

    if data[line_counter][5]:
        ax.annotate("", xy=(x + 2, y + .5), xytext=(x + 1.25, y + .5),
                    arrowprops=dict(width=3, fc='black')
                    )
        probability = 100 if data[line_counter][-1] is None else data[line_counter][7]
        ax.text(x + 1.625, y + .25,
                str(probability) + '%', fontname='Cambria',
                fontsize=15, color="black",
                ha="center", va="center"
                )
        ax.text(x + 1.625, y + .7,
                ("{:.1e}".format(data[line_counter][0]) + ' s'),
                fontsize=15, color="black",
                ha="center", va="center", fontname='Cambria'
                )
        ax.text(x + 1.625, y + .95, str(data[line_counter][4]) + ' decay',
                fontsize=15, color="black", ha="center",
                va="center", style='italic', fontname='Cambria'
                )
        chain_img(data, line_counter + 1, x + 2.25, y, ax, depth + 1)
    else:
        ax.text(x + .5, y + .25, data[line_counter][4], fontsize=13,
                color="black", ha="center", va="center", fontname='Cambria'
                )
        plt.plot(x, y)


def activity_calc(decay_const, time, A0):
    """Activity of every chain member at ``time``, given ``A0`` of the parent.

    Solves dN/dt = A N for the bidiagonal decay matrix A by matrix exponential,
    rather than by the Bateman closed form.

    The closed form divides by the products of (lambda_j - lambda_i) for
    j != i, so it fails outright the moment two chain members share a decay
    constant -- which happens whenever two members have the same tabulated
    half-life, and database half-lives are rounded. Nudging the duplicates
    apart does not rescue it: the terms then scale as 1/epsilon and
    1/epsilon**2, and with three equal constants that needs ~1e18 of dynamic
    range, so double precision cancels away every significant figure and the
    answer can even come out negative.

    expm handles repeated eigenvalues exactly, and agrees with the closed form
    to ~1e-12 when the decay constants are distinct.
    """
    lam = np.asarray(decay_const, dtype=float)
    size = len(lam)

    matrix = np.zeros((size, size))
    np.fill_diagonal(matrix, -lam)          # each member decays away
    if size > 1:
        matrix[1:, :-1] += np.diag(lam[:-1])  # ...and feeds the next one

    populations = np.zeros(size)
    populations[0] = A0 / lam[0]
    return list(lam * (expm(matrix * time) @ populations))


def activity_img(data, A0, results_dir):
    """Activity of each chain member against time, on log-log axes."""
    decay_const = []
    titles = []
    summ = 0
    t_arr = []
    for line in data:
        if line[4] == 'STABLE':
            break
        t_arr = np.concatenate([t_arr,
                                np.linspace(
                                            summ, summ + line[0] * 13,
                                            700, endpoint=False
                                            )])
        summ += line[0] * 10
        decay_const.append(np.log(2) / line[0])
        titles.append(line[3])

    if not titles:
        return None

    act = [[] for _ in range(len(titles))]
    for t in t_arr:
        for i, A in enumerate(activity_calc(decay_const, t, A0)):
            act[i].append(A)

    use_style('classic')
    figure, ax = plt.subplots(figsize=(11.07, 8), facecolor='white')
    for i, isotope in enumerate(act):
        plt.plot(t_arr, isotope, label=titles[i], linewidth=2)

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Time, s', fontsize=15, fontname='Cambria')
    plt.ylabel('Activity, Bq', fontsize=15, fontname='Cambria')
    plt.xticks(fontsize=12, fontname='Cambria')
    plt.yticks(fontsize=12, fontname='Cambria')
    plt.ylim(bottom=A0 / 1E8, top=A0 * 10)
    for line in data:
        if line[0] and line[0] > 1e4:
            plt.xlim(left=1)
            break

    # prop must be a FontProperties/dict; a bare string is not accepted.
    plt.legend(prop={'family': 'Cambria', 'size': 22})
    ax.grid(visible=True)
    ax.set_facecolor('white')

    path = results_dir / 'activity.png'
    plt.savefig(path, bbox_inches='tight')
    plt.close(figure)
    return path


# --------------------------------------------------------------------------- #
# application
# --------------------------------------------------------------------------- #
def load_chart(cursor):
    """Ground-state half-life and decay type of every nuclide in the database."""
    cursor.execute('''SELECT DecayModes.hl_sec, DecayModes.isotope_z,
                   DecayModes.isotope_n, Isotope.name, DecayTypes.name
                   FROM DecayModes join Isotope join DecayTypes
                   where DecayModes.isotope_z = Isotope.z
                   and DecayModes.isotope_n = Isotope.n
                   and DecayTypes.id = DecayModes.decaytype_id
                   and DecayModes.e_level_mev=0
                   order by DecayModes.isotope_z asc,
                   DecayModes.isotope_n asc,
                   DecayModes.probability asc;'''
                   )
    return cursor.fetchall()


def start(user_input, frame1, frame2, frame3, notebook, cursor, config):
    """Rebuild all three figures for the isotope named in the entry box."""
    for frame in (frame2, frame3):
        for child in frame.winfo_children():
            child.destroy()

    data = load_chart(cursor)

    xy = None
    for line in data:
        if user_input == line[3]:
            xy = [line[2], line[1]]
            break
    if xy is None:
        messagebox.showerror('Unknown isotope',
                             f'{user_input!r} is not in the database.\n'
                             'Enter a name such as 81Ga.')
        return

    results_dir = HERE / config['paths']['results_dir']
    results_dir.mkdir(parents=True, exist_ok=True)

    heatmap_path = heatmap_img(data, xy, results_dir)

    chain = []
    isotope_chain_data(xy[1], xy[0], chain, cursor)
    if not chain:
        messagebox.showerror('No decay data',
                             f'No decay modes recorded for {user_input}.')
        return

    activity_path = activity_img(chain, config['view'].getfloat('initial_activity'),
                                 results_dir)

    use_style('classic')
    fig, ax = plt.subplots(figsize=(len(chain) * 3, 7), facecolor='white')
    ax.set_aspect('equal', 'box')
    ax.set_facecolor("whitesmoke")
    ax.axis('off')
    chain_img(chain, 0, 1, 1, ax)
    chain_path = results_dir / 'chain.png'
    plt.savefig(chain_path, bbox_inches='tight')
    plt.close(fig)

    notebook.tab(frame2, state='normal')
    notebook.tab(frame3, state='normal')

    heatmap_pic = Image.open(heatmap_path).resize((960, 640),
                                                  Image.Resampling.LANCZOS)
    chain_pic = Image.open(chain_path)
    if chain_pic.size[0] > 960:
        scale = int(chain_pic.size[1] * 960 / chain_pic.size[0])
        chain_pic = chain_pic.resize((960, scale), Image.Resampling.LANCZOS)

    heatmap_pic = ImageTk.PhotoImage(heatmap_pic)
    chain_pic = ImageTk.PhotoImage(chain_pic)

    label_heatmap = tk.Label(frame1, image=heatmap_pic)
    label_chain = tk.Label(frame2, image=chain_pic)
    # Tk does not keep a Python reference to a PhotoImage, so without these the
    # images are garbage-collected when start() returns and the tabs go blank.
    label_heatmap.image = heatmap_pic
    label_chain.image = chain_pic

    label_heatmap.place(relx=0, rely=0.995, anchor='sw')
    label_chain.place(relx=.5, rely=.5, anchor='center')

    if activity_path is not None:
        activity_pic = ImageTk.PhotoImage(
            Image.open(activity_path).resize((960, 640),
                                             Image.Resampling.LANCZOS))
        label_activity = tk.Label(frame3, image=activity_pic)
        label_activity.image = activity_pic
        label_activity.place(relx=.5, rely=.5, anchor='center')


def main():
    config = load_config()
    password = (os.environ.get('ISOTOPEDB_GUEST_PASSWORD')
                or getpass.getpass('Password for the read-only account: '))

    nudat = mysql.connector.connect(
        host=config['database']['host'],
        user=config['database']['user'],
        password=password,
        database=config['database']['database'],
        )
    cursor = nudat.cursor()

    try:
        root = tk.Tk()
        root.geometry("960x720")

        notebook = ttk.Notebook(root)
        notebook.pack()

        frm_map = tk.Frame(notebook, height=720, width=960, bg='white')
        frm_chain = tk.Frame(notebook, height=720, width=960, bg='white')
        frm_activity = tk.Frame(notebook, height=720, width=960, bg='white')

        for frame in (frm_map, frm_chain, frm_activity):
            frame.pack(fill='both', expand=1)

        notebook.add(frm_map, text='Isotope map')
        notebook.add(frm_chain, text='Decay chains', state='disabled')
        notebook.add(frm_activity, text='Activity plot', state='disabled')

        label1 = tk.Label(frm_map, text='Isotope:', bg='white', font='Cambria')
        label1.place(relx=0.05, rely=0.04, anchor='center')

        entry1 = tk.Entry(frm_map, width=6, bg='whitesmoke')
        entry1.insert(0, config['view']['default_isotope'])
        entry1.place(relx=0.12, rely=0.04, anchor='center')

        btn = tk.Button(frm_map, text='confirm',
                        command=(lambda: start(entry1.get(), frm_map,
                                               frm_chain, frm_activity,
                                               notebook, cursor, config
                                               ))
                        )
        btn.place(relx=0.17, rely=0.04, anchor='center')

        root.mainloop()
    finally:
        cursor.close()
        nudat.close()


if __name__ == '__main__':
    main()
