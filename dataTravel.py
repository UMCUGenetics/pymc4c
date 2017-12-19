import sys
import numpy as np
import collections

import matplotlib.pyplot as pp
import bisect as bs

class DataTraveller:

	selection = collections.defaultdict(list)
	def __init__(self, filePath):
		reStOut = np.load(filePath)
		try:
			self.byRead = reStOut['byread'].item()
		except KeyError:
			self.byRead = reStOut['byreads'].item()
		self.byRegion = reStOut['byregion'].item()

		self.history = []

		tmp = dict()
		for key in self.byRegion:

			sub = dict()
			for i,reg in enumerate(self.byRegion[key]):
				#print reg
				sub[i] = (reg[0],[x[0] for x in reg[1]])
				# print i,reg
				# ~exit()
			tmp[key] = sub#range(len(self.byRegion[key]))
		self.curRegion = tmp
		print 'curRegion',self.curRegion['chr8'][0]

		self.curRead = self.byRead.keys()

	def coordsToRegionIndex(self, chromosome, start, end):
		regionList = self.curRegion[chromosome]
		print 'curRegion',self.curRegion['chr8']
		#print chromosome, start, end
		#print regionList
		# Use bisect implementation to quickly find a matching position
		left = bs.bisect([regionList[x][0][1] for x in regionList],start)
		right = bs.bisect([regionList[x][0][0] for x in regionList],end)
		#print left,right
		refLen = len(regionList)-1

		# Don't bother beyond the last position in the list
		right = min(right, refLen)
		return left, right


	# [(chr,regionID),...]
	def getReadsByRegions(self, regions):
		targetReads = set()
		for regionIndex in regions:
			region = self.byRegion[regionIndex[0]][regionIndex[1]]
			for read in region[1]:
				targetReads.add(read)
		return targetReads


	# [readID,...]
	# Returns dict[chr][]
	def getRegionsByReads(self, reads):
		selection = collections.defaultdict(list)
		for read in reads:
			for location in self.byRead[read]:
				if location[0] not in selection:
					selection[location[0]] = collections.defaultdict(list)
				selection[location[0]][location[1]].append(read) #byread[read]
		return selection


	def filterByRegions(self, regions):
		mergeSet = set()
		for region in regions:
			# print region
			chrom = region[0]
			location = region[1]
			mergeSet = mergeSet.union(self.byRegion[chrom][location])
		return mergeSet


	def applyAnchor(self, regions):
		self.history.append(regions)
		
		#print self.curRegion['chr8']
		print self.getReadsByRegions(regions)
		allReads = [x[0] for x in self.getReadsByRegions(regions)]
		
		self.curRead = [x for x in self.curRead if x in allReads]
		self.curRegion = self.getRegionsByReads(self.curRead)
		
		#print self.curRegion['chr8']

	def getSelectedInfo(self,chromosomes=None):
		infoList = []

		if chromosomes == None:
			chromosomes = self.curRegion.keys()

		for chromosome in chromosomes:
			# print chromosome
			for region in self.curRegion[chromosome]:
				# print region
				#infoList.append(byRegion[chromosome][region])
				#print 'derp',region
				#print 'Do iets',region
				infoList.append(
					(
					region,
					chromosome,
					self.byRegion[chromosome][region][0][0],
					self.byRegion[chromosome][region][0][1],
					len(self.byRegion[chromosome][region][1]),
					len([x for x in self.byRegion[chromosome][region][1] if x[0] in self.curRead]),
					)
				)
				
		infoList.sort()

		#return [(x[0],x[1][0],x[1][1],x[2],x[3]) for x in infoList]

		#return [(x[0]+':'+str(x[1][0])+'-'+str(x[1][1]),x[2]) for x in infoList]
		return infoList


	def getAppliedAnchors(self):
		anchors = collections.defaultdict(set)
		for group in self.history:
			for region in group:
				anchors[region[0]].add(region[1])

		anchorCounts = dict()
		for chromosome in self.byRegion.keys():
			if chromosome in anchors:
				anchorCounts[chromosome] = len(anchors[chromosome])
			else:
				anchorCounts[chromosome] = 0

		return anchorCounts


	def getChromosomeRegions(self):
		matches= dict()
		for chromosome in self.byRegion.keys():
			if chromosome in self.curRegion:
				matches[chromosome] = len(self.curRegion[chromosome])
			else:
				matches[chromosome] = 0

		return matches
