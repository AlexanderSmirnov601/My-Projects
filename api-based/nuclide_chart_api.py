"""Nuclide chart viewer backed by the IAEA Live Chart API.

The same three-tab Tkinter viewer as ``nuclide-chart-viewer`` -- half-life
heatmap with inset, decay chain, chain activity -- but reading nuclide data
from the public IAEA nuclear data service instead of the IsotopeDB MySQL
database. The figure code is imported from the original viewer, so both
versions render identically; only the data layer differs.

No database or credentials are needed. The one bulk API request is cached on
disk, so the viewer also works offline after the first run.
"""

import importlib.util
import sys
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'Results'
DEFAULT_ISOTOPE = '81Ga'
INITIAL_ACTIVITY_BQ = 1E+9
MIN_BRANCH_PERCENT = 1.0        # second decay branches below this are ignored

sys.path.insert(0, str(HERE))
import iaea_client

# Change in (Z, N) produced by each IAEA decay-mode label. Modes with no
# single child nuclide (fission) map to None.
MODE_CHANGES = {
    'B-': (1, -1), 'B+': (-1, 1), 'EC': (-1, 1), 'EC+B+': (-1, 1),
    'A': (-2, -2), 'IT': (0, 0),
    'B-N': (1, -2), 'B-2N': (1, -3), 'B-A': (-1, -3),
    'ECP': (-2, 1), 'B+P': (-2, 1),
    'N': (0, -1), '2N': (0, -2), 'P': (-1, 0), '2P': (-2, 0),
    '2B-': (2, -2), '2EC': (-2, 2), '2B+': (-2, 2),
    'SF': None, 'ECSF': None,
}


def load_viewer():
    """Import the database viewer module and reuse its figure functions."""
    path = HERE.parent / 'nuclide-chart-viewer' / 'nuclide_chart_viewer.py'
    spec = importlib.util.spec_from_file_location('nuclide_chart_viewer', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules['nuclide_chart_viewer'] = module
    spec.loader.exec_module(module)
    return module


class NuclideData:
    """Chart and decay-chain data derived from the IAEA ground-states table."""

    def __init__(self):
        states = iaea_client.ground_states()
        # Branching-ratio columns hold '%', which itertuples() cannot expose.
        states = states.rename(columns={'decay_1_%': 'decay_1_pct',
                                        'decay_2_%': 'decay_2_pct'})
        states['name'] = ((states.z + states.n).astype(int).astype(str)
                          + states.symbol.astype(str))
        self.states = states
        self.by_zn = {(int(row.z), int(row.n)): row
                      for row in states.itertuples()}

    def chart_rows(self):
        """(hl_sec, z, n, name, decay_type) per nuclide, like load_chart()."""
        rows = []
        for row in self.states.itertuples():
            stable = str(row.half_life) == 'STABLE'
            if not stable and pd.isna(row.half_life_sec):
                continue            # no measured half-life
            decay = 'STABLE' if stable else str(row.decay_1)
            hl = None if stable else float(row.half_life_sec)
            rows.append((hl, int(row.z), int(row.n), row.name, decay))
        return rows

    def _decay_row(self, state, mode, percent):
        """One chain row shaped like the database query result."""
        z, n = int(state.z), int(state.n)
        change = MODE_CHANGES.get(mode)
        child_z = child_n = None
        if change is not None:
            child_z, child_n = z + change[0], n + change[1]
        probability = None if pd.isna(percent) else float(percent)
        return (float(state.half_life_sec), z, n, state.name, mode,
                child_z, child_n, probability, 0.0)

    def chain(self, z, n, data=None, visited=None):
        """Walk the decay chain from (z, n), main branch plus one side branch."""
        data = [] if data is None else data
        visited = set() if visited is None else visited
        if (z, n) in visited or len(visited) >= 200:
            return data
        visited.add((z, n))

        state = self.by_zn.get((z, n))
        if state is None:
            return data
        if str(state.half_life) == 'STABLE' or pd.isna(state.half_life_sec):
            data.append((None, z, n, state.name, 'STABLE',
                         None, None, None, 0.0))
            return data

        main = self._decay_row(state, str(state.decay_1),
                               state.decay_1_pct)
        data.append(main)
        if main[5] is not None:
            self.chain(main[5], main[6], data, visited)
        else:
            return data

        second_mode = str(state.decay_2)
        second_percent = state.decay_2_pct
        if (second_mode not in ('nan', '')
                and not pd.isna(second_percent)
                and float(second_percent) >= MIN_BRANCH_PERCENT
                and MODE_CHANGES.get(second_mode) is not None
                and MODE_CHANGES[second_mode] != (0, 0)):
            branch = self._decay_row(state, second_mode, second_percent)
            data.append(branch)
            self.chain(branch[5], branch[6], data, visited)
        return data


def start(user_input, frame1, frame2, frame3, notebook, viewer, data):
    """Rebuild all three figures for the isotope named in the entry box."""
    for frame in (frame2, frame3):
        for child in frame.winfo_children():
            child.destroy()

    chart = data.chart_rows()
    xy = None
    for line in chart:
        if user_input == line[3]:
            xy = [line[2], line[1]]
            break
    if xy is None:
        messagebox.showerror('Unknown isotope',
                             f'{user_input!r} is not in the IAEA table.\n'
                             'Enter a name such as 81Ga.')
        return

    RESULTS.mkdir(exist_ok=True)
    heatmap_path = viewer.heatmap_img(chart, xy, RESULTS)

    chain = data.chain(xy[1], xy[0])
    if not chain:
        messagebox.showerror('No decay data',
                             f'No decay modes recorded for {user_input}.')
        return

    activity_path = viewer.activity_img(chain, INITIAL_ACTIVITY_BQ, RESULTS)

    viewer.use_style('classic')
    fig, ax = plt.subplots(figsize=(len(chain) * 3, 7), facecolor='white')
    ax.set_aspect('equal', 'box')
    ax.set_facecolor('whitesmoke')
    ax.axis('off')
    viewer.chain_img(chain, 0, 1, 1, ax)
    chain_path = RESULTS / 'chain.png'
    plt.savefig(chain_path, bbox_inches='tight')
    plt.close(fig)

    notebook.tab(frame2, state='normal')
    notebook.tab(frame3, state='normal')

    from PIL import Image, ImageTk
    heatmap_pic = ImageTk.PhotoImage(
        Image.open(heatmap_path).resize((960, 640), Image.Resampling.LANCZOS))
    chain_pic = Image.open(chain_path)
    if chain_pic.size[0] > 960:
        scale = int(chain_pic.size[1] * 960 / chain_pic.size[0])
        chain_pic = chain_pic.resize((960, scale), Image.Resampling.LANCZOS)
    chain_pic = ImageTk.PhotoImage(chain_pic)

    label_heatmap = tk.Label(frame1, image=heatmap_pic)
    label_chain = tk.Label(frame2, image=chain_pic)
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
    viewer = load_viewer()
    data = NuclideData()

    root = tk.Tk()
    root.geometry('960x720')
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
    entry1.insert(0, DEFAULT_ISOTOPE)
    entry1.place(relx=0.12, rely=0.04, anchor='center')
    btn = tk.Button(frm_map, text='confirm',
                    command=(lambda: start(entry1.get(), frm_map, frm_chain,
                                           frm_activity, notebook, viewer,
                                           data)))
    btn.place(relx=0.17, rely=0.04, anchor='center')

    root.mainloop()


if __name__ == '__main__':
    main()
