'''global properties of how to construct a graph'''

Path2ViterbiJar = '../Script/ViterbiJar/ViterbiAlgorithm.jar'
CutOffContact = 7.5
G_ConstructType = 2
NumberOfClusterInClusteringGraph = 20
'''
Construct type:
    @0: DisMx(pair) <= threshold, single edges, weights as number of contacts
    @1: min(pair) <= threshold, single edge, has weight as exp(-var)
    @2: max(pair) <= threshold, single edge, has weight as exp(-var)
'''

List_Colors = ['darkgreen', 'darkred', 'darkblue', 'lightblue','lightgreen', 'lightcoral','yellow', 'purple','green', 'red', 'cyan', 'magenta',
               'pink', 'brown', 'orange', 'teal', 'coral', 'blue', 'lime',
               'lavender', 'turquoise', 'tan', 'salmon', 'gold',
               'lightpurple', 'black']