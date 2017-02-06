### --- Picking up at S06 RFS01 here --- ###
# TODO: Check correctness of this implementation, might contain some basepair shifts due to 0vs1 based genome stuff

import m2p
import prep
from Bio.Seq import Seq
import numpy as np
import pandas as pd

import pysam
import parse
import collections

fastaIdFormat='>{3};pt{0}/{2}\n{1}\n'
referenceNameFormat='RD:{};IN:{}'

def splitStringTo(item,maxLen=50):
	return [str(item[ind:ind+maxLen]) for ind in range(0, len(item), maxLen/2)]


# In case of leftover sequence the last bit is attached to the preceeding chunk
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
			i+1,
			val,
			len(split),
			baseId)
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
	#item = df.loc[args.id]

	leftSeq, rightSeq = getPrimerSeqs(dataInfo.loc[args.id])
	writePrimerFasta(leftSeq, rightSeq, args.output)

	#for item in df.id:
	#	print 'Working on',item
		#leftSeq, rightSeq = getPrimerSeqs(df.loc[item])



# This class can perhaps be replaced by a panda frame later
class SimpleRead(object):
	#fastaIdFormat='PR{}_{}-{}'#s06.fastaIdFormat
	def __init__(self, read, prmLen):
		# May want to replace this bit in the future: 'PR{}_{}-{}' -> fastaIdFormat
		#primerType,part,parts=list(parse.parse('PR{}_{}-{}',read.query_name))
		primerType,part,parts=list(parse.parse('PR{};pt{}/{}',read.query_name))

		self.prmType = primerType
		self.prmFlag = read.flag & ~(1<<8) # Set 9th bit to 0 to state primary alignment
		#self.readID = parse.parse('SN{}_{}',read.reference_name)[0]
		self.readID = parse.parse(referenceNameFormat,read.reference_name)[0]
		self.startAln = read.reference_start
		self.endAln = read.reference_start + read.infer_query_length(always=True)
		self.prmSize = prmLen[int(primerType)-1] # Should be length of the original primer?
						# prm_len{prm_info(ai,1)}(prm_info(ai,8))
		self.alnErr = [sum(x[1] for x in read.cigartuples if x[0] in [1,2])]
		self.winInd = [int(part)]
		self.winNum = int(parts)
		self.cigar = [read.cigarstring]

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

# This stuff happens after the botwie part in s06


def combinePrimers(insam,prmLen,qualThreshold=.20):
	samfile = pysam.AlignmentFile(insam, "rb")

	mySet = set()
	myDict = collections.defaultdict(list)
	for read in samfile:
		myDict[parse.parse(referenceNameFormat,read.reference_name)[0]].append(SimpleRead(read,prmLen))

		# Just to speed up testing this implementation:
		#if len(myDict) > 5000:
		#	break

	#samfile.reset()
	sortedKeys = myDict.keys()
	sortedKeys.sort(key=int)
	simpleList = []

	def checkPrimers(matchList):
		matchList.sort(key=lambda x: x.startAln)
		primTypeDict = collections.defaultdict(int)
		for x in matchList:
			primTypeDict[x.prmType] += 1 #toList()

		cutList = []

		if matchList == []:
			return cutList

		#if max(primTypeDict[x] for x in primTypeDict)>1: # Otherwise no problems found
		print ''

		if matchList[0].prmFlag&(1<<4)==16:
			print 'Accept: First bit'
			cutList.append([0,matchList[0].endAln])
		for i in xrange(0,len(matchList)-1):
			#print matchList[i].toList()
			if matchList[i].prmFlag&(1<<4)==0 and \
				matchList[i+1].prmFlag&(1<<4)==16 and \
				matchList[i].prmType != matchList[i+1].prmType:
				print 'Accept:',i,i+1
				cutList.append([matchList[i].startAln,matchList[i+1].endAln])

		#print matchList[-1].toList()
		if matchList[-1].prmFlag&(1<<4)==0:
			print 'Accept: Last bit'
			cutList.append([matchList[-1].startAln,None])

		for i,x in enumerate(matchList):
			print i,x.toList()

		for x in cutList:
			print x

		return cutList

	for key in sortedKeys:
		# Remove low quality matches
		myDict[key] =  [x for x in myDict[key] if (x.alnErr[0] / float(x.endAln-x.startAln) < qualThreshold)]

		#if len(myDict[key]) < 5:
		#	continue
		for side in set(x.prmType for x in myDict[key]):
			thisSide = [x for x in myDict[key] if x.prmType.startswith(side)]
			if len(thisSide) > 1:
				# Merge latter listed reads into earlier reads if overlapping
				for j in xrange(len(thisSide)-1, -1, -1):
					for i in xrange(j-1, -1, -1):
						curRead = thisSide[i]
						nextRead = thisSide[j]

						if curRead.endAln + 10 >= nextRead.startAln \
								and curRead.startAln <= nextRead.endAln + 10 \
								and curRead.prmFlag&(1<<4) == nextRead.prmFlag&(1<<4): # Are on same strand
							#print side,i,j,":",curRead.startAln,curRead.endAln,"---",nextRead.startAln,nextRead.endAln
							curRead.startAln = min(curRead.startAln,nextRead.startAln)
							curRead.endAln = max(curRead.endAln,nextRead.endAln) #nextRead.endAln
							curRead.winInd.extend(nextRead.winInd)
							curRead.alnErr.extend(nextRead.alnErr)
							curRead.cigar.extend(nextRead.cigar)
							myDict[key].remove(nextRead)
							break

		checkPrimers(myDict[key])

		for x in myDict[key]:
			#print x.tostring()
			simpleList.append(x.toList())

	pdfColumns = ['prmType','prmFlag','readID',
		'startAln','endAln','prmSize',
		'alnErr','winInd','winNum','cigar']
	#print len(simpleList),len(simpleList[0])
	pdFrame = pd.DataFrame(simpleList,columns=pdfColumns)

	return pdFrame


def cleaveReads(args):
	dataInfo = m2p.load(args.info_file)
	primerLen1 = len(dataInfo['pr1_seq'][args.id])
	primerLen2 = len(dataInfo['pr2_seq'][args.id])
	print [primerLen1,primerLen2]
	prmInfo = combinePrimers(args.rfs_file,[primerLen1,primerLen2])
	exit()
	print prmInfo
	#print prmInfo['prmType']
	# TODO: This appears to be a switch that should be an argument upon calling the script
	if dataInfo['seq_plt'][args.id] == 'PacBio':
		lowQualIndex = prmInfo['alnErr']/prmInfo['prmSize'] > 0.15
	elif dataInfo['seq_plt'][args.id] == 'NanoPore':
		lowQualIndex = prmInfo['alnErr']/prmInfo['prmSize'] > 0.20
	else:
		print 'Unknown platform:',dataInfo['seq_plt'][args.id]
		exit(1)

	print lowQualIndex
	#print len(lowQualIndex),sum(lowQualIndex)
	#print prmInfo
