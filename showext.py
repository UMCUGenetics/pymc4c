import sys
import numpy as np
import collections

import matplotlib.pyplot as pp

selection = collections.defaultdict(list)

restout = np.load(sys.argv[1])
byread = restout['byread'].item()
byregion = restout['byregion'].item()

targetRegions = set()
targetReads = set()

for key in byregion:
    for region in byregion[key]:
        if len(region[1])>50:
            # print key,region[0]
            for read in region[1]:
                targetReads.add(read)
            targetRegions.add((key,region[0]))

for read in targetReads:
    for location in byread[read]:
        if location[0] not in selection:
            selection[location[0]] = collections.defaultdict(list)
        selection[location[0]][location[1]].append(read) #byread[read]

for reg in targetRegions:
    print reg
print '\n'
# print selection

values = []
labels = []

sortedKeys = selection.keys()
sortedKeys.sort()
for chromosome in sortedKeys:
    sortedRegions = selection[chromosome].keys()
    sortedRegions.sort()
    for region in sortedRegions:
        reads = selection[chromosome][region]
        if len(reads) > 25:
            values.append(len(reads))
            print region,chromosome,byregion[chromosome][region][0],len(reads),(chromosome,byregion[chromosome][region][0]) in targetRegions
        #if len(region)>3:
        #    print region,len(selection[region])


pp.bar(range(len(values)),values,log=True)
#pp.xticks(np.arange(len(values)),labels)

pp.show()
