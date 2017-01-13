import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot(inf, outf, samples, legend=True):

    df = pd.read_table(inf, index_col=0, skiprows=1).transpose()
    df = df[df.index.isin(samples)]
    df = df.reindex([samples])
    labels  = samples

    columns = df.columns

    sorted_columns = []
    for c in columns:
        sorted_columns.append([df[c].sum() ,c])
    sorted_columns = sorted(sorted_columns)
    columns = [c[1] for c in sorted_columns]

    N = len(samples)
    colors = plt.cm.Set3(np.linspace(0, 1, len(columns)))
    colors = colors[::-1]
    ind = np.arange(N)
    width = 1
    fig, ax = plt.subplots()
    counts = []
    for smpl in samples:
        counts.append(df.at[smpl, columns[0]])

    bars = ax.bar(ind, counts, width, color=colors[0])
    b = counts
    for br in bars:
        l = columns[0].split(';')[-1]
        if l == '':
            l = 'NA'
        br.set_label(l)

    i = 1
    for col in columns[1:]:
        counts = []
        for smpl in samples:
            counts.append(df.at[smpl, col])
        bars = ax.bar(ind, counts, width, bottom=b, color=colors[i])
        for br in bars:
            l = col.split(';')[-2]+' '+col.split(';')[-1]
            if l == '':
                l = col.split(';')[-2]
            if l == '':
                l = col.split(';')[-3]
            if l == '':
                l='NA'
            br.set_label(l)
        b = [b[j]+counts[j] for j in range(N)]
        i+=1
    plt.xticks(ind+width/2., labels,  rotation=90)
    plt.ylim([0,1])
    plt.xlim([0,N+6])
    ax.tick_params(labelsize=32)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles[::(-1*N)][:20], labels[::(-1*N)][:20], fontsize=20)

    # Write the figure
    fig = plt.gcf()
    fig.set_size_inches(30,20)
    fig.savefig(outf, bbox_inches='tight')
    #plt.close('all')

colors = {name: plt.get_cmap(name) for name in [
    'Blues', 'BuGn', 'BuPu','GnBu', 'Greens',
    'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn',
    'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn',
    'YlGnBu', 'YlOrBr', 'YlOrRd'
]}

def plot2(inf, outf, colorkey='Phylum', colormap=None, numcatcol=6, width=30,height=20,tickfont=22):
    df = pd.read_table(inf)

    df = df.sort_values(colorkey)

    groupsizes = {}

    for x, row in df.iterrows():
        try:
            groupsizes[row[colorkey]] += 1
        except:
            groupsizes[row[colorkey]] = 1

    colorlist = []

    for group in groupsizes:
        colorlist += list(colors[colormap[group]](np.linspace(0.1,0.9,groupsizes[group])))

    fig, ax = plt.subplots()

    N = len(df.columns[numcatcol:])
    ind = np.arange(N)
    width = 1

    b = None
    first = True
    for x, row in df.iterrows():
        if first:
            bars = ax.bar(ind, row.tolist()[numcatcol:], width, color=colorlist.pop(0))
            b = row.tolist()[numcatcol:]
            first=False
            continue
        bars = ax.bar(ind, row.tolist()[numcatcol:], width, bottom=b, color=colorlist.pop(0))
        for i in range(len(b)):
            b[i] += row.tolist()[numcatcol:][i]

    plt.xticks(ind+width/2., df.columns[numcatcol:],  rotation=90)
    plt.ylim([0,1])
    plt.xlim([0,N+6])
    ax.tick_params(labelsize=tickfont)


    # Write the figure
    fig = plt.gcf()
    fig.set_size_inches(width,height)
    fig.savefig(outf, bbox_inches='tight')