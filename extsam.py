import sys
import pysam
import numpy as np
import bisect as bs

def mapToRefSite(refSiteList,mappedPos):
    # Use bisect implementation to quickly find a matching position
    pos = bs.bisect_left(refSiteList,mappedPos)
    refLen = len(refSiteList)-1

    # Don't bother beyond the last position in the list
    pos = min(pos, refLen)

    # Move start and end positions to match our heuristic
    left = pos
    right = pos
    while refSiteList[left][0] >= mappedPos[0] + 10 and left > 0:
        left -= 1
    while refSiteList[right][1] <= mappedPos[1] - 10 and right < refLen:
        right += 1

    return [left, right]

# reflist = [[100,104],[110,120],[140,144],
#             [200,204],[210,220],[240,244],
#             [300,304],[310,320],[340,344]]
# placelist = [[130,200],[205,255],[295,360],[270,280]]
# for thisplace in placelist:
#     #thisplace=[205,255]
#     print mapToRefSite(reflist,thisplace)
#exit()
print 'Loading restrsites, this takes a while...'
restrefs=np.load(sys.argv[2])['restrsites'].item()
print 'Finished loading, moving on'

insam = sys.argv[1]
samfile = pysam.AlignmentFile(insam, "rb")

prevRead = samfile.next()
prevResult = [-1,-1]
curStack = []

for read in samfile:
    if not read.is_unmapped:
        if read.reference_name not in restrefs:
            continue
        result = mapToRefSite(restrefs[read.reference_name],[read.reference_start, read.reference_start + read.infer_query_length(always=True)])

        # If two subsequent reads:
            # were mapped to the same chromosome,
            # and to the same strand,
            # and with at least one overlapping restriction site
            # and have the same mother read...
        if prevRead.reference_id == read.reference_id \
            and prevRead.is_reverse == read.is_reverse \
            and result[0] <= prevResult[1] and result[1] >= prevResult[0] \
            and int(dict(item.split(":") for item in read.query_name.split(";"))['RD']) \
                == int(dict(item.split(":") for item in prevRead.query_name.split(";"))['RD']):
                 curStack.append((read,result))
        else:
            if curStack != []:
                for i in range(min([x[1][0] for x in curStack]), max([x[1][1] for x in curStack])+1):
                    readID = int(dict(item.split(":") for item in prevRead.query_name.split(";"))['RD']),[x[1] for x in curStack]
                    restrefs[read.reference_name][i].append(readID)
                #print int(dict(item.split(":") for item in prevRead.query_name.split(";"))['RD']),[x[1] for x in curStack], min([x[1][0] for x in curStack]), max([x[1][1] for x in curStack])
            curStack = [(read,result)]

        prevRead = read
        prevResult = result

for key in restrefs:
    curList = restrefs[key]
    keyLen = len(curList)


    # Insert one item to remove some edge issues
    curList.insert(0,[0,0])

    # Turn restriction site locations into meaningful regions
    for i in xrange(0,keyLen):
        curList[i][0] = sum(curList[i][:2])/2
        curList[i][1] = sum(curList[i+1][:2])/2

    # Remove empty values from the dataset
    for x in xrange(keyLen,0,-1):
        if len(curList[x]) <= 2:
            curList.pop(x)
        #if max(len(curList[x]),len(curList[x+1])) == 2:
        #    curList.pop(x+1)
    # Remove the added edge item
    #curList.pop(0)
    print keyLen, len(curList)
    for x in curList:
        print x[:2],len(x)-2
