### --- Picking up at S06 RFS01 here --- ###
# TODO: Check correctness of this implementation, might contain some basepair shifts due to 0vs1 based genome stuff

import m2p
import prep
from Bio.Seq import Seq
import numpy as np
import pandas as pd

import re
import pysam
import parse
import collections

fastaIdFormat='>PI:{};WI:{};WN:{}\n{}\n'
referenceNameFormat='>RD:{};IN:{}'

def splitStringTo(item,maxLen=50):
	return [str(item[ind:ind+maxLen]) for ind in range(0, len(item), maxLen/2)]


# In case of leftover sequence the last bit is attached to the preceding chunk
def seqToFasta(sequence,baseId):
	targetLen = 50
	split = splitStringTo(sequence,maxLen=targetLen)
	if len(split) > 1:
	 	while len(split[0]) != len(split[-1]):
			popper = split.pop(-1)
		split[-1] = split[-1][:targetLen/2]+popper
	outString = ''
	for i,val in enumerate(split):
		outString+=fastaIdFormat.format(
			baseId,
			i+1,
			len(split),
			val
			)
	return outString


# TODO: Improve error message on faulty primes sequences
def getPrimerSeqs(dataInfo):
	leftSeq = prep.getFastaSequence(
			dataInfo.genome_build,
			dataInfo.vp_chr,
			dataInfo.pr1_pos[0]-300,
			dataInfo.pr1_pos[1]).upper()
	leftIndex = leftSeq.rfind(dataInfo.re1_seq)
	leftPrimerSeq = Seq(leftSeq[leftIndex:]).reverse_complement()
	assert leftPrimerSeq.find(dataInfo.pr1_seq) >= 0, 'Primer sequence is wrong\n'+str(dataInfo)

	rightSeq = prep.getFastaSequence(
			dataInfo.genome_build,
			dataInfo.vp_chr,
			dataInfo.pr2_pos[0]-1,
			dataInfo.pr2_pos[1]+300).upper()
	rightIndex = rightSeq.find(dataInfo.re1_seq) + len(dataInfo.re1_seq) - 1
	rightPrimerSeq = rightSeq[:rightIndex];
	assert rightPrimerSeq.find(dataInfo.pr2_seq) >= 0, 'Primer sequence is wrong\n'+str(dataInfo)

	return leftPrimerSeq,rightPrimerSeq


def writePrimerFasta(leftPrimerSeq,rightPrimerSeq,targetFile):
	with open(targetFile,'w') as outFasta:
		outFasta.write(seqToFasta(leftPrimerSeq,'PR1'))
		outFasta.write(seqToFasta(rightPrimerSeq,'PR2'))


def makePrimerFasta(args):
	dataInfo = m2p.load(args.infile)

	leftSeq, rightSeq = getPrimerSeqs(dataInfo.loc[args.id])
	writePrimerFasta(leftSeq, rightSeq, args.output)


# This class can perhaps be replaced by a panda frame later
class SimpleRead(object):
	#fastaIdFormat='PR{}_{}-{}'#s06.fastaIdFormat
	def __init__(self, read, prmLen):
		primerDict = dict(item.split(":") for item in read.query_name.split(";"))

		self.prmType = primerDict['PI'][-1]
		self.prmFlag = read.flag & ~(1<<8) # Set 9th bit to 0 to state primary alignment
		self.readID=int(dict(item.split(":") for item in read.reference_name.split(";"))['RD'])
		self.startAln = read.reference_start
		self.endAln = read.reference_start + read.infer_query_length(always=True)
		self.prmSize = prmLen[int(self.prmType)-1] # Should be length of the original primer?
						# prm_len{prm_info(ai,1)}(prm_info(ai,8))
		self.alnErr = [sum(x[1] for x in read.cigartuples if x[0] in [1,2])]
		self.winInd = [int(primerDict['WI'])]
		self.winNum = int(primerDict['WN'])
		self.cigar = [read.cigarstring]

		# Cigar string information
		# M 	BAM_CMATCH 		0
		# I 	BAM_CINS 		1
		# D 	BAM_CDEL 		2
		# N 	BAM_CREF_SKIP 	3
		# S 	BAM_CSOFT_CLIP 	4
		# H 	BAM_CHARD_CLIP 	5
		# P 	BAM_CPAD 		6
		# = 	BAM_CEQUAL 		7
		# X 	BAM_CDIFF 		8

	def toList(self):
		return [
			self.prmType,
			self.prmFlag,
			self.readID,
			self.startAln,
			self.endAln,
			self.prmSize,
			self.alnErr,
			self.winInd,
			self.winNum,
			self.cigar]

	def tostring(self):
		return '\t'.join(str(x) for x in [
			self.prmType,
			self.prmFlag,
			self.readID,
			self.startAln,
			self.endAln,
			self.prmSize,
			self.alnErr,
			self.winInd,
			self.winNum,
			self.cigar])

# This stuff happens after the bowtie2 part in s06
def groupPrimers(matchList):
	# Split by primer type to ensure we only combine primes of the same sort
	for side in set(x.prmType for x in matchList):
		thisSide = [x for x in matchList if x.prmType.startswith(side)]
		if len(thisSide) > 1:
			# Merge latter listed reads into earlier reads if overlapping
			for j in xrange(len(thisSide)-1, -1, -1):
				for i in xrange(j-1, -1, -1):
					curRead = thisSide[i]
					nextRead = thisSide[j]

					# Allow a slight offset to consider overlapping to identify
					# win_i and win_i+2 correctly, and check if both are on the
					# same strand
					if curRead.endAln + 10 >= nextRead.startAln \
							and curRead.startAln <= nextRead.endAln + 10 \
							and curRead.prmFlag&(1<<4) == nextRead.prmFlag&(1<<4):
						# Merge information
						curRead.startAln = min(curRead.startAln,nextRead.startAln)
						curRead.endAln = max(curRead.endAln,nextRead.endAln)
						curRead.winInd.extend(nextRead.winInd)
						curRead.alnErr.extend(nextRead.alnErr)
						curRead.cigar.extend(nextRead.cigar)
						# Remove last item in list of merged information parts
						matchList.remove(nextRead)
						# Merge further in subsequent loops, not in this one
						break

	return matchList


def findCuts(matchList):
	# Ensure the list provided is sorted by start positions
	matchList.sort(key=lambda x: x.startAln)
	cutList = []

	if matchList == []:
		return cutList

	# If the first primer points to the start, accept it
	if matchList[0].prmFlag&(1<<4)==16:
		cutList.append([None,matchList[0].endAln])

	# Any two subsequent primers pointing toward eachother are accepted,
	# but ensure they are not the same primer
	for i in xrange(0,len(matchList)-1):
		if matchList[i].prmFlag&(1<<4)==0 and \
				matchList[i+1].prmFlag&(1<<4)==16 and \
				matchList[i].prmType != matchList[i+1].prmType:
			cutList.append([matchList[i].startAln,matchList[i+1].endAln])
		# TODO: Check if we need to add +1 to the end position above
		# If we want to remove reads where 2 of the same primer type are
		# pointing toward eachother we could return an empty list here

	# If the last primer points to the end, accept it
	if matchList[-1].prmFlag&(1<<4)==0:
		cutList.append([matchList[-1].startAln,None])

	return cutList


def combinePrimers(insam,prmLen,qualThreshold=.20):
	samfile = pysam.AlignmentFile(insam, "rb")

	# Create a dict with lists of SimpleRead using referenceName as key
	myDict = collections.defaultdict(list)
	for read in samfile:
		myDict[int(dict(item.split(":") for item in read.reference_name.split(";"))['RD'])].append(SimpleRead(read,prmLen))

	# Put data in referenceName sorted order
	sortedKeys = myDict.keys()
	sortedKeys.sort(key=int)
	cutAllList = []

	for key in sortedKeys:
		# Remove low quality matches
		matchedRead = [x for x in myDict[key] if (x.alnErr[0] / float(x.endAln-x.startAln) < qualThreshold)]

		# Group overlapping primers on a read
		grouped = groupPrimers(matchedRead)

		# Determine where cuts should be based on mapped primers
		cutAllList.append([key,findCuts(grouped)])

	return cutAllList


def applyCuts(inFile,outFile,cutList,cutDesc='PC'):
	readId=-1
	readName=''
	readSeq=''
	with open(inFile,'r') as faFile, open(outFile,'w') as dumpFile:
		# Indexing is lead by cutlist, contains less than or equal to inFile
		for cut in cutList:
			# No cuts can be made, ignore this sequence
			if cut[1] == []:
				continue

			# Assuming both lists are sorted by read id (int),
		 	# play catchup between the two lists
			while cut[0] > readId:
				readName=faFile.next().rstrip()
				readSeq=faFile.next().rstrip()
				readId=int(dict(item.split(":") for item in readName[1:].split(";"))['RD'])

			# Both lists are now either aligned or inFile is ahead
			if readId == cut[0]:
				# Split sequence, dump information
				for i,x in enumerate(cut[1]):
					# Extend the identifier
					dumpFile.write(readName+';'+cutDesc+':'+str(i)+'\n')#+':'+str(x[0])+'-'+str(x[1])
					# Dump the actual sub sequence
					dumpFile.write(readSeq[x[0]:x[1]]+'\n')


def cleaveReads(args):
	dataInfo = m2p.load(args.info_file)
	primerLen1 = len(dataInfo['pr1_seq'][args.id])
	primerLen2 = len(dataInfo['pr2_seq'][args.id])
	prmCuts = combinePrimers(args.rfs_file,[primerLen1,primerLen2])
	applyCuts(args.infasta,args.outfasta,prmCuts)


def findRestrictionSeqs(inFile,outFile,restSeqs,cutDesc='RC'):
	compRestSeqs = [str(Seq(x).reverse_complement()) for x in restSeqs]
	restSeqs.extend(compRestSeqs)
	reSeqs='|'.join(restSeqs)
	cutList = []

	with open(inFile,'r') as faFile, open(outFile,'w') as dumpFile:
		for read in faFile:
			readName=read.rstrip() #faFile.next().rstrip()
			readSeq=faFile.next().rstrip()

			matches = [[x.start(), x.end()] for x in (re.finditer(reSeqs, readSeq))]
			thisCut = []
			if matches != []:
				thisCut.append([None,matches[0][1]])
				for i in xrange(len(matches)-1):
					thisCut.append([matches[i][0],matches[i+1][1]])
				thisCut.append([matches[-1][0],None])

			# Split sequence, dump information
			for i,x in enumerate(thisCut):
				# Extend the identifier
				dumpFile.write(readName+';'+cutDesc+':'+str(i)+'\n')#+':'+str(x[0])+'-'+str(x[1])
				# Dump the actual sub sequence
				dumpFile.write(readSeq[x[0]:x[1]]+'\n')

			cutList.append([readName,thisCut])
	return cutList


def splitReads(args):
	dataInfo = m2p.load(args.info_file)
 	restSeqs = [dataInfo['re1_seq'][args.id],dataInfo['re2_seq'][args.id]]
	# TODO: Substitute reference genome with reads (?)
	findRestrictionSeqs(args.infasta,args.outfasta,restSeqs)
