import sys
import re
from Bio.Seq import Seq
import numpy as np


def findSites(refFile,restSeqs,lineLen=50):
    compRestSeqs = [str(Seq(x).reverse_complement()) for x in restSeqs]
    restSeqs.extend(compRestSeqs)
    reSeqs='|'.join(restSeqs)
    print reSeqs
    restSitesDict = dict()

    with open(refFile,'r') as reference:
        curChrom = None
        offset = -lineLen
        matches = []
        readSeq = ''
        for line in reference:
            if line[0] == '>':
                matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start()])
                readSeq = 'N'*lineLen*2
                offset = -lineLen
                matches = []
                restSitesDict[line[1:].rsplit()[0]] = matches
            else:
                readSeq = readSeq[lineLen:]+line.rsplit()[0].upper()
                matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start() < lineLen])
                offset += lineLen
        matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start()])

    return restSitesDict

np.savez_compressed(sys.argv[2],restrsites=findSites(sys.argv[1],['GACC']))
