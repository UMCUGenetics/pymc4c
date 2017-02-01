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

def splitStringTo(item,maxLen=50):
	return [str(item[ind:ind+maxLen]) for ind in range(0, len(item), maxLen/2)]


# In case of leftover sequence the last bit is attached to the preceeding chunk
def seqToFasta(sequence,baseId):
	split = splitStringTo(sequence)
	if len(split) > 1 and len(split[0]) != len(split[-1]):
		split[len(split)-2] += split.pop(-1)
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
		primerType,part,parts=list(parse.parse('PR{}_{}-{}',read.query_name))

		self.prmType = primerType
		self.prmFlag = read.flag & ~(1<<8) # Set 9th bit to 0 to state primary alignment
		self.readID = parse.parse('SN{}_{}',read.reference_name)[0]
		self.startAln = read.reference_start
		self.endAln = read.reference_start + read.infer_query_length(always=True)
		self.prmSize = prmLen[int(primerType)-1] # Should be length of the original primer?
						# prm_len{prm_info(ai,1)}(prm_info(ai,8))
		self.alnErr = sum(x[1] for x in read.cigartuples if x[0] in [1,2])
		self.winInd = [part]
		self.winNum = parts
		self.cigar = [read.cigarstring]

		# M 	BAM_CMATCH 	0
		# I 	BAM_CINS 	1
		# D 	BAM_CDEL 	2
		# N 	BAM_CREF_SKIP 	3
		# S 	BAM_CSOFT_CLIP 	4
		# H 	BAM_CHARD_CLIP 	5
		# P 	BAM_CPAD 	6
		# = 	BAM_CEQUAL 	7
		# X 	BAM_CDIFF 	8

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


def combinePrimers(insam,prmLen):
	samfile = pysam.AlignmentFile(insam, "rb")

	mySet = set()
	myDict = collections.defaultdict(list)
	for read in samfile:
		myDict[parse.parse('SN{}_{}',read.reference_name)[0]].append(SimpleRead(read,prmLen))

		# Just to speed up testing this implementation:
		if len(myDict) > 5000:
			break

	#samfile.reset()
	sortedKeys = myDict.keys()
	sortedKeys.sort(key=int)
	simpleList = []

	for key in sortedKeys:
		for side in set(x.prmType for x in myDict[key]):
			thisSide = [x for x in myDict[key] if x.prmType.startswith(side)]
			if len(thisSide) > 1:
				for i in xrange(len(thisSide)-2, -1, -1):
					curRead = thisSide[i]
					nextRead = thisSide[i+1]

					# Should we allow a little mismatch here?
					# It does influence results, check isOLEx() in matlab?
					if curRead.endAln >= nextRead.startAln:
						curRead.endAln = nextRead.endAln
						curRead.winInd.extend(nextRead.winInd)
						curRead.alnErr += nextRead.alnErr
						curRead.cigar.extend(nextRead.cigar)
						myDict[key].remove(nextRead)

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
	prmInfo = combinePrimers(args.rfs_file,[primerLen1,primerLen2])

	print prmInfo['prmType']
	# TODO: This appears to be a switch that should be an argument upon calling the script
	if dataInfo['seq_plt'][args.id] == 'PacBio':
		lowQualIndex = prmInfo['alnErr']/prmInfo['prmSize'] > 0.15
	elif dataInfo['seq_plt'][args.id] == 'NanoPore':
		lowQualIndex = prmInfo['alnErr']/prmInfo['prmSize'] > 0.20
	else:
		print 'Unknown platform:',dataInfo['seq_plt'][args.id]
		exit(1)

	print sum(lowQualIndex)
