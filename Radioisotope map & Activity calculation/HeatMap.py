import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import tkinter as tk
from tkinter import ttk
from PIL import ImageTk, Image
import mysql.connector
import getpass


def start(user_input, frame1, frame2, frame3, notebook, cursor):
    for child in frame2.winfo_children():
        child.destroy()
        for child in frame3.winfo_children():
            child.destroy()

    cursor.execute('''use IsotopeDB;'''
                   )
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
    data = cursor.fetchall()

    xax, yax = 180, 120
    arr = [[None for x in range(xax)] for y in range(yax)]
    annotation = [[None for x in range(xax)] for y in range(yax)]

    for row in data:
        if row[0] is None:
            arr[row[1]][row[2]] = 1E+30
            annotation[row[1]][row[2]] = row[3] + '\n' + str(row[4])
        else:
            arr[row[1]][row[2]] = row[0]
            annotation[row[1]][row[2]] = row[3] + '\n' + str(row[4])
            annotation[row[1]][row[2]] += '\n' + "{:.2e}".format(row[0]) + ' s'

    for x in range(len(arr)):
        for y in range(len(arr[x])):
            if arr[x][y] is None:
                arr[x][y] = np.nan

    log_norm = LogNorm(vmin=1E-3, vmax=1E+15)

    ax, fig = plt.subplots(figsize=[21, 15])
    plt.style.use('seaborn')

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

    plt.xlim([1, 178])
    plt.ylim([0, 118])
    for line in data:
        if user_input == line[3]:
            xy = [line[2], line[1]]
            break

    axins = ax.inset_axes([0.595, 0.058, 0.45, 0.53])
    for z, line in enumerate(annotation):
        for n, atom in enumerate(line):
            if z < xy[1]-2 or z > xy[1]+2 or n < xy[0]-2 or n > xy[0]+2:
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

    if xy[0] >= 2:
        axins.set_xlim(xy[0]-2, xy[0]+3)
    else:
        axins.set_xlim(0, xy[0]+3)

    if xy[1] >= 2:
        axins.set_ylim(xy[1]-2, xy[1]+3)
    else:
        axins.set_ylim(0, xy[1]+3)

    axins.set_xticklabels(axins.get_xmajorticklabels(), fontsize=14,
                          weight='bold', fontname='Cambria'
                          )
    axins.set_yticklabels(axins.get_ymajorticklabels(), fontsize=14,
                          weight='bold', fontname='Cambria'
                          )

    patch, pp1, pp2 = mark_inset(ax, axins,
                                 loc1=1, loc2=1, edgecolor='firebrick',
                                 linewidth=3, alpha=1
                                 )
    pp1.loc1 = 2
    pp2.loc1 = 2
    if xy[0] < 120:
        pp1.loc2 = 1
        pp2.loc2 = 1
    else:
        pp1.loc2 = 3
        pp2.loc2 = 3

    plt.savefig('Results/HeatMap.png',
                bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    data = []
    isotope_chain_data(xy[1], xy[0], data, cursor)
    set_length = len(data)

    plt.style.use('classic')
    activity_img(data, 1E+9)

    fig, ax = plt.subplots(figsize=(set_length*3, 7), facecolor='white')
    ax.set_aspect('equal', 'box')
    ax.set_facecolor("whitesmoke")
    ax.axis('off')
    chain_img(data, 0, 1, 1, ax)
    plt.savefig('Results/chain.png',
                bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    notebook.tab(frm_chain, state='normal')
    notebook.tab(frm_activity, state='normal')

    heatmap_pic = Image.open('Results/HeatMap.png')
    chain_pic = Image.open('Results/chain.png')
    activity_pic = Image.open('Results/activity.png')

    heatmap_pic = heatmap_pic.resize((960, 640), Image.Resampling.LANCZOS)
    if chain_pic.size[0] > 960:
        scale = int(chain_pic.size[1]*960/chain_pic.size[0])
        chain_pic = chain_pic.resize((960, scale),
                                     Image.Resampling.LANCZOS
                                     )
    activity_pic = activity_pic.resize((960, 640), Image.Resampling.LANCZOS)

    heatmap_pic = ImageTk.PhotoImage(heatmap_pic)
    chain_pic = ImageTk.PhotoImage(chain_pic)
    activity_pic = ImageTk.PhotoImage(activity_pic)

    label_heatmap = tk.Label(frame1, image=heatmap_pic)
    label_chain = tk.Label(frame2, image=chain_pic)
    label_activity = tk.Label(frame3, image=activity_pic)

    label_heatmap.pack()
    label_chain.pack()
    label_activity.pack()

    label_heatmap.place(relx=0, rely=0.995, anchor='sw')
    label_chain.place(relx=.5, rely=.5, anchor='center')
    label_chain.place(relx=.5, rely=.5, anchor='center')

    nudat.commit()
    return


def chain_img(data, line_counter, x, y, ax):
    ax.add_patch(Rectangle((x, y), 1, 1, fill=True,
                           ec='black', fc='skyblue', linewidth=2
                           )
                 )
    ax.text(x+.5, y+.5, data[line_counter][3],
            fontsize=15, fontname='Cambria',
            color="black", ha="center",
            va="center", weight='bold'
            )
    for line in data[line_counter+1:]:
        if (
            line[3] == data[line_counter][3]
            and data[data.index(line)-1][4] == 'STABLE'
        ):
            ax.annotate("", xy=(x+.5, y-1), xytext=(x+.5, y-.25),
                        arrowprops=dict(width=3, fc='black'),
                        fontname='Cambria'
                        )
            ax.text(x+.25, y-.625,
                    str(data[data.index(line)][7])+'%',
                    fontsize=15, color="black", ha="center",
                    va="center", rotation=-90, fontname='Cambria'
                    )
            ax.text(x+.7, y-.625,
                    "{:.1e}".format(data[data.index(line)][0])+' s',
                    fontsize=15, color="black", ha="center",
                    va="center", rotation=-90, fontname='Cambria'
                    )
            ax.text(x+.95, y-.625, str(data[data.index(line)][4])+' decay',
                    color="black", ha="center", va="center", style='italic',
                    rotation=-90, fontsize=15, fontname='Cambria'
                    )
            plt.plot(x, y-2.25)
            chain_img(data, data.index(line)+1, x, y - 2.25, ax)

            break

    if data[line_counter][5]:
        ax.annotate("", xy=(x+2, y+.5), xytext=(x+1.25, y+.5),
                    arrowprops=dict(width=3, fc='black')
                    )
        if data[line_counter][-1] is None:
            none = 100
        else:
            none = data[line_counter][7]
        ax.text(x+1.625, y+.25,
                str(none)+'%', fontname='Cambria',
                fontsize=15, color="black",
                ha="center", va="center"
                )
        ax.text(x+1.625, y+.7,
                ("{:.1e}".format(data[line_counter][0]) + ' s'),
                fontsize=15, color="black",
                ha="center", va="center", fontname='Cambria'
                )
        ax.text(x+1.625, y+.95, str(data[line_counter][4])+' decay',
                fontsize=15, color="black", ha="center",
                va="center", style='italic', fontname='Cambria'
                )
        x += 2.25
        line_counter += 1
        chain_img(data, line_counter, x, y, ax)
    else:
        ax.text(x+.5, y+.25, data[line_counter][4], fontsize=13,
                color="black", ha="center", va="center", fontname='Cambria'
                )
        plt.plot(x, y)
        return


def activity_img(data, A0):
    decay_const = []
    titles = []
    summ = 0
    t_arr = []
    for line in data:
        if line[4] == 'STABLE':
            break
        else:
            t_arr = np.concatenate([t_arr,
                                    np.linspace(
                                                summ, summ+line[0]*13,
                                                700, endpoint=False
                                                )])
            summ += line[0]*10
            decay_const.append(np.log(2)/line[0])
            titles.append(line[3])

    # t_arr = np.linspace(0, time, 3000, endpoint=True)
    act = [[] for i in range(len(titles))]

    for t in t_arr:
        calc = activity_calc(decay_const, t, A0)
        for i, A in enumerate(calc):
            act[i].append(A)

    figure, ax = plt.subplots(figsize=(11.07, 8), facecolor='white')
    for i, isotope in enumerate(act):
        plt.plot(t_arr, isotope, label=titles[i], linewidth=2)

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Time, s', fontsize=15, fontname='Cambria')
    plt.ylabel('Activity, Bq', fontsize=15, fontname='Cambria')
    plt.xticks(fontsize=12, fontname='Cambria')
    plt.yticks(fontsize=12, fontname='Cambria')
    plt.ylim(bottom=A0/1E8)
    plt.ylim(top=A0*10)
    for line in data:
        if line[0] and line[0] > 1e4:
            plt.xlim(left=1)
            break

    plt.legend(prop='Cambria', fontsize=22)
    ax.grid(visible=True)
    ax.set_facecolor('white')
    plt.savefig('Results/activity.png',
                bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def activity_calc(decay_const, time, A0):
    A = []
    eq1 = np.arange(len(decay_const))
    eq2 = np.arange(len(decay_const))
    for k, h in enumerate(decay_const):
        eq2 = []
        for i, h1 in enumerate(decay_const[:k+1]):
            eq1 = []
            for j, h2 in enumerate(decay_const[:k+1]):
                eq1.append(decay_const[j]-h1)
            eq1.pop(eq1.index(0))
            eq2.append(np.exp(-h1*time)/np.prod(eq1))
        act = A0*np.prod(decay_const[:k+1])/decay_const[0]*np.sum(eq2)
        A.append(act)
    return A


def isotope_chain_data(z, n, data, cursor):
    cursor.execute(
                    '''use IsotopeDB;'''
                    )
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
    nudat.commit()
    if arr:
        if arr[0][5]:
            data.append(arr[0])
            isotope_chain_data(arr[0][5], arr[0][6], data, cursor)
        elif arr[0][4] == 'STABLE' and len(arr) >= 2 and arr[1][4] == 'IT':
            data.append(arr[1])
            data.append(arr[0])
            return data
        else:
            data.append(arr[0])
            return data
        if len(arr) >= 2 and arr[1][7] > 1 and arr[1][-1] == 0:
            data.append(arr[1])
            isotope_chain_data(arr[1][5], arr[1][6], data, cursor)
    else:
        return data


nudat = mysql.connector.connect(
    host='34.67.94.150',
    user='guest',
    password = 'password',
    )
cursor = nudat.cursor()

root = tk.Tk()
root.geometry("960x720")

notebook = ttk.Notebook(root)
notebook.pack()

frm_map = tk.Frame(notebook, height=720, width=960, bg='white')
frm_chain = tk.Frame(notebook, height=720, width=960, bg='white')
frm_activity = tk.Frame(notebook, height=720, width=960, bg='white')

frm_map.pack(fill='both', expand=1)
frm_chain.pack(fill='both', expand=1)
frm_activity.pack(fill='both', expand=1)

notebook.add(frm_map, text='Isotope map')
notebook.add(frm_chain, text='Decay chains', state='disabled')
notebook.add(frm_activity, text='Activity plot', state='disabled')


label1 = tk.Label(frm_map, text='Isotope:', bg='white', font='Cambria')
label1.pack()
label1.place(relx=0.05, rely=0.04, anchor='center')

entry1 = tk.Entry(frm_map, width=6, bg='whitesmoke')
entry1.pack()
entry1.insert(0, '81Ga')
entry1.place(relx=0.12, rely=0.04, anchor='center')

btn = tk.Button(frm_map, text='confirm',
                command=(lambda: start(entry1.get(), frm_map,
                                       frm_chain, frm_activity,
                                       notebook, cursor
                                       ))
                )
btn.place(relx=0.17, rely=0.04, anchor='center')

root.mainloop()
cursor.close()
nudat.close()