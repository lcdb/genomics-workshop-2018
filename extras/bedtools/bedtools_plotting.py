import pybedtools
from pybedtools.contrib import plotting

import time
import os
import pybedtools
from pybedtools.contrib import plotting
from matplotlib import pyplot as plt
import seaborn as sns

outdir = 'images'

if not os.path.exists(outdir):
    os.makedirs(outdir)

colors = ['k', '#004573', '#cc0000']

def fix_args(x):
    if isinstance(x, str):
        x = x.split()
    g = []
    skipnext = False
    for i,j  in zip(x[:-1], x[1:]):
        if i.startswith('-') and not j.startswith('-'):
            j = os.path.basename(j)
            g.append(' '.join([i, j]))
            skipnext = True
        else:
            if not skipnext:
                g.append(os.path.basename(i))
            skipnext = False

        def sort(x):
            if x == 'bedtools':
                return 0
            if x.startswith('-a'):
                return 2
            if x.startswith('-b'):
                return 3
            if x.startswith('-'):
                return 9
            return 1

    return ' '.join(sorted(g, key=sort))

def plot_bedtool(filea, fileb=None, method=None, **kwargs):
    """
    Use for BEDTools programs that use -a and -b input arguments.  Filenames
    `a` and `b` are used for `method` of the BedTool class, and `kwargs` are
    sent to that method.

    The result is a plot of `a`, `b`, and the result, with the commandline
    argument as the plot title.
    """
    a = pybedtools.BedTool(filea)
    if fileb is not None:
        b = pybedtools.BedTool(fileb)
        kwargs['b'] = fileb
        result = getattr(a, method)(**kwargs)
        bfn = b.fn
        bts = [result, b, a]
        labels = ['result', os.path.basename(b.fn), os.path.basename(a.fn)]
    else:
        result = getattr(a, method)(**kwargs)
        bfn = None
        bts = [result, a]
        labels = ['result', os.path.basename(a.fn)]

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    ybase = 0
    yheight = 1
    ylabels = []
    yticks = []
    pad = 0.4
    any_stranded = sum(sum(j.strand != '.' for j in i) for i in bts)
    for color, bt, label in zip(colors, bts, labels):
        if bt is None:
            continue
        ax.axhline(ybase, color='k', zorder=100)
        ybase += pad

        if any_stranded:
            strands = ['+', '-', '.']
            strand_labels = {
                '+': '(+)',
                '-': '(-)',
                '.': '( )',
            }
        else:
            strands = ['.']
            strand_labels = {'.': ''}
        for strand in strands:
            sbt = bt.filter(lambda x: x.strand == strand)
            ylabels.append(label + ' ' + strand_labels[strand])
            track = plotting.Track(
                    sbt, visibility='squish', alpha=0.5, ybase=ybase, color=color, stranded=False)
            yticks.append(track.midpoint)
            ybase = track.ymax + pad
            ax.axhline(ybase, color='0.7', linewidth=0.3)
            if strand != '.':
                ybase += pad
            ax.add_collection(track)
    ax.axhline(ybase, color='k', zorder=100)
    ax.set_yticklabels(ylabels, fontdict=dict(name='Monospace'))
    ax.set_yticks(yticks)

    cmds = result._cmds
    cmds[0] = 'bedtools ' + pybedtools.settings._prog_names[cmds[0]]
    cleaned_cmds = fix_args(cmds)
    ax.set_title(cleaned_cmds, fontdict=dict(name='Monospace'))
    ax.axis('tight')
    fig.subplots_adjust(top=0.8, bottom=0.15)
    sns.despine(ax=ax, left=True, right=True, top=True)
    plt.tick_params(axis='y', which='both', left='off')
    ax.axis(ymin=0)
    return ax, cleaned_cmds

a = pybedtools.BedTool('x.bed')
b = pybedtools.BedTool('y.bed')

def s(x):
    x.strand = '.'
    return x

with open('genome.chromsizes', 'w') as fout:
    fout.write('chr1\t5000\n')


def plot_and_save(*args, **kwargs):
    ax, cmds = plot_bedtool(*args, **kwargs)
    fig = ax.figure
    fn = cmds.replace(' ', '_') + '.png'
    fig.savefig(os.path.join(outdir, fn))

plot_and_save(filea=a, fileb=b, method='intersect')
plot_and_save(filea=a, fileb=b, u=True, method='intersect')
plot_and_save(filea=b, fileb=a, u=True, method='intersect')
plot_and_save(filea=a, fileb=b, v=True, method='intersect')
plot_and_save(filea=b, fileb=a, v=True, method='intersect')
plot_and_save(filea=a, method='merge')
plot_and_save(filea=a, fileb=b, method='subtract')
plot_and_save(filea=a, method='flank', l=10, r=0, g='genome.chromsizes')
plot_and_save(filea=a, method='slop', b=50, g='genome.chromsizes')



plt.show()
