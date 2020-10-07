def WriteNumericalList2File(Path2OutFile, Lines):
    thefile = open(Path2OutFile, 'w')
    for line in Lines:
        for l in line:
            thefile.write("%s\t" % l)
        thefile.write("\n")
    thefile.close()
def WriteList2File(Path2OutFile, Lines):
    thefile = open(Path2OutFile, 'w')
    for line in Lines:
        thefile.write("%s\n" %line)
    thefile.close()
def WriteMap2File(Path2OutFile, Maps):
    Lines = []
    for K in Maps.keys():
        line = str(K) + "\t" + str(Maps[K])
        Lines.append(line)
    WriteList2File(Path2OutFile, Lines)
def ReadClusterIds(Path2File):
    f = open ( Path2File , 'r')
    MyMap = {}
    for line in f:
        arr = map(int,line.split())
        MyMap.update({arr[0]:arr[1:]})
    Arr = []
    for i in range(len(MyMap)):
        Arr.append(MyMap[i])
    return Arr
def ReadGraphFeatureFile(Path2File):
    f = open ( Path2File , 'r')
    Map_VertexID_Score = {}
    Map_EdgeID_Score = {}
    for line in f:
        if line.startswith("#"):
            continue
        arr = line.split()
        if len(arr)==3:
            Map_VertexID_Score.update({int(arr[0]):[float(arr[1]), float(arr[2])]})
        elif len(arr)==6:
            Map_EdgeID_Score.update({(int(arr[0]), int(arr[1])):[float(arr[2]),float(arr[3]),float(arr[4]),float(arr[5])]})
    return (Map_VertexID_Score, Map_EdgeID_Score)
def ReadLines(Path2File):
    f = open ( Path2File , 'r')
    Lines = []
    for line in f:
        if len(line) > 0:
            Lines.append(line)
    return Lines

def ReadViterbiOutFile(Path2File):
    from collections import OrderedDict
    f = open(Path2File, 'r')
    MyMap = {}
    for line in f:
        arr = list(map(int, line.split()))
        MyMap.update({arr[0]: arr[1]})
    OD = OrderedDict(sorted(MyMap.items()))
    return OD.values()



        