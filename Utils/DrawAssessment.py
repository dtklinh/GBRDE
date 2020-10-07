from plot_bars import plot_domains
import numpy as np, seaborn as sns
from matplotlib.pyplot import figure
from csb.io import load
#sns.set(style='ticks', palette='Set1', context='notebook', font_scale=1.)
sns.set(palette='Set2', context='notebook', font_scale=1.5)



def load_membership_linh(serial):
    rs = open('../data/Membership/'+str(serial)+'.txt').readlines()[0]
    mbs = rs.split('\t')
    mbs = np.array([int(c) for c in mbs])
    return mbs

# if False:
#     #proteins = ('pyruvate', 'T7', 'GroEL', 'AST', 'AdH')
#     proteins = ('T7',)
#     #serial = [2733, 441, 1537, 344, 243]
#     serial = [441]
#     for i, name in enumerate(proteins):
#         new_ = []
#         membership = load('../data/{}.pkl'.format(name))
#         membership1 = load_membership_linh(serial[i])
#         new_.append(membership1)
#         new_.append(membership[0])
#         for j in membership[3:]:
#             new_.append(j)
#             print len(membership1), len(j)
#
#         from csb.io import dump
#
#         dump(new_, '../data/linh/{}.pkl'.format(name))

    #load Linh results


if True:
    order = (
    'Literature', 'Ours', 'Spectrus', 'Spectrus (K=2)', 'Spectrus (K=3) *', 'DynDom')

    method_dict = {'AdH': ['Ours', 'DynDom', 'Spectrus (K=2)', 'Spectrus (K=3) *'],
                   'AST': ['Ours', 'DynDom', 'Spectrus'],
                   'pyruvate': ['Ours', 'DynDom', 'Spectrus', 'Literature'],
                   'GroEL': ['Ours', 'DynDom', 'Spectrus', 'Literature'],
                   'T7': [ 'Ours','Literature', 'Spectrus', 'DynDom']}

    names = {'T7': 'T7 RNA polymerase',
             'pyruvate': 'Pyruvate phosphate dikinase',
             'GroEL': 'GroEL',
             'AdH': 'Alcohol dehydrogenase',
             'AST': 'Aspartate aminotransferase'}

    colors = sns.color_palette()


    def sort_membership(membership):

        new_membership = np.zeros(len(membership), 'i')
        labels = list(set(membership.tolist()))
        counts = map(membership.tolist().count, labels)
        labels = np.array(labels)[np.argsort(counts)[::-1]]

        for i, label in enumerate(labels):
            new_membership += i * (membership == label).astype('i')

        return new_membership


    def add_plot(ax, name, ref_name='Literature'):

        membership = map(sort_membership, load('../data/linh/{}.pkl'.format(name)))
        methods = method_dict[name]
        indices = map(order.index, methods)
        indices = np.argsort(indices)[::-1]
        membership = map(membership.__getitem__, indices)
        methods = map(methods.__getitem__, indices)

        K = max(map(max, membership))
        if ref_name in methods:
            plot_domains(K, membership, colors, methods, 0.55, methods.index(ref_name), ax=ax)
        else:
            plot_domains(K, membership, colors, methods, 0.55, ax=ax)
        ax.annotate(names[name], xy=(0., 0.65), fontsize=16, xycoords='axes fraction')
        ax.grid(False)
        ax.set_ylim(-0.5, 5.5)

    fig = figure(figsize=(8, 12))

    proteins = ('pyruvate', 'T7', 'GroEL', 'AST', 'AdH')
    n = len(proteins)

    ref_names = {'HIV': 'Gibbs (prior 2)'}

    for counter, protein in enumerate(proteins, 1):
        ax = fig.add_subplot(n, 1, counter)
        # ax.set_facecolor('Black')

        add_plot(ax, protein, ref_name=ref_names.get(protein, 'Literature'))

    fig.savefig('../latex/img/fig11.pdf', bbox_inches='tight')