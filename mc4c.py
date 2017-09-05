import argparse
import sys

import log
import numpy as np
import mc4ctools as mc
import pandas as pd


def makePrimerFasta(args):
	settings = mc.loadIni(args.inifile)
	primerSeqs = mc.getPrimerSeqs(settings)
	mc.writePrimerFasta(primerSeqs, args.outfile)


def cleaveReads(args):
	settings = mc.loadIni(args.inifile)
	primerLens = [len(x) for x in settings['prm_seq']]
	primers = ['']
	primers.extend(settings['prm_seq'])
	#print primers
	prmCuts = mc.combinePrimers(args.bamfile,primerLens)
	#print prmCuts[:10]
	mc.applyCuts(args.fastqfile,args.outfile,prmCuts,primers)


def splitReads(args):
	settings = mc.loadIni(args.inifile)
	restSeqs = settings['re_seq']
	# TODO: Substitute reference genome with reads (?)
	mc.findRestrictionSeqs(args.fastqfile,args.outfile,restSeqs)


def findRefRestSites(args):
	settings = mc.loadIni(args.inifile)
	restSeqs = settings['re_seq']
	restDict = mc.findReferenceRestSites(args.fastafile,restSeqs)
	np.savez_compressed(args.restfile,restrsites=restDict)


def exportToPlot(args):
	settings = mc.loadIni(args.inifile)
	print 'Loading restrsites, this takes a while...'
	restrefs=np.load(args.restfile)['restrsites'].item()
	print 'Finished loading, moving on'
	byRegion,byRead,pdFrame = mc.exportToPlot(settings,restrefs,args.bamfile)

	#dupSet = mc.findDuplicates(settings,byRead,byRegion)
	#pdFrame['Duplicate'] = np.where(pdFrame['CircleId'].isin(dupSet), True, False)

	#print pdFrame
	np.savez_compressed(args.plotfile,
		pdframe=pdFrame,
		pdcolumns=pdFrame.columns,
		pdindex=pdFrame.index)

	np.savez_compressed(args.plotfile+'_extra',
		byregion=byRegion,
		byread=byRead)


def markDuplicates(args):
	settings = mc.loadIni(args.inifile)
	exFile = np.load(args.extra)

	try:
		byRead = exFile['byread'].item()
	except KeyError:
		byRead = exFile['byreads'].item()
	byRegion = exFile['byregion'].item()

	pdFile = np.load(args.pdframe)
	pdFrame = pd.DataFrame(pdFile['pdframe'],columns=pdFile['pdcolumns'],index=pdFile['pdindex'])
	dupSet = mc.findDuplicates(settings,byRead,byRegion)

	#df['dup']=np.where(pd.Series(df.index).isin([1,5]),True,False)
	#pdFrame['Duplicate'] = np.where(pdFrame['CircleId'].isin(dupSet), True, False)

	pdFrame['Duplicate'] = np.where(pd.Series(pdFrame.index).isin(dupSet), True, False)

	np.savez_compressed(args.outfile,
		pdframe=pdFrame,
		pdcolumns=pdFrame.columns,
		pdindex=pdFrame.index)

# Huge wall of argparse text starts here
def main():
	descIniFile = 'File containing experiment specific details'
	descFqFile = 'Fastq file containing actual data from sequencing'

	parser = argparse.ArgumentParser(
		description="MC4C pipeline for processing multi-contact data")
	subparsers = parser.add_subparsers()

	#
	parser_mkprfa = subparsers.add_parser('makeprimerfa',
		description='Make a fasta file of primer sequences')
	parser_mkprfa.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_mkprfa.add_argument('outfile',
		type=str,
		help='Fasta file with primer sequences')
	parser_mkprfa.set_defaults(func=makePrimerFasta)

	#
	parser_clvprm = subparsers.add_parser('cleavereads',
		description='Cleave reads by primer sequences')
	parser_clvprm.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_clvprm.add_argument('bamfile',
		type=str,
		help='Bam file that is created after makeprimerfa results were mapped by bowtie2')
	parser_clvprm.add_argument('fastqfile',
		type=str,
		help=descFqFile)
	parser_clvprm.add_argument('outfile',
		type=str,
		help='Fastq file to dump primer cleaved sequences into')
	parser_clvprm.set_defaults(func=cleaveReads)

	#
	parser_splrest = subparsers.add_parser('splitreads',
		description='Split reads by restriction site sequences')
	parser_splrest.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_splrest.add_argument('fastqfile',
		type=str,
		help=descFqFile)
	parser_splrest.add_argument('outfile',
		type=str,
		help='Fasta file to dump restriction site split sequences into')
	parser_splrest.set_defaults(func=splitReads)

	#
	parser_refrest = subparsers.add_parser('refrestr',
		description='Determine restriction site coordinates on reference by their sequences')
	parser_refrest.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_refrest.add_argument('fastafile',
		type=str,
		help='Fasta file of reference genome containing all chromosomes (eg hg19)')
	parser_refrest.add_argument('restfile',
		type=str,
		help='Numpy compressed file containing restriction site coordinates')
	parser_refrest.set_defaults(func=findRefRestSites)

	#
	parser_export = subparsers.add_parser('export',
		description='Combine and export results for interactive plotting')
	parser_export.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_export.add_argument('bamfile',
		type=str,
		help='Bam file after mapping previously split fastq data to reference genome (eg hg19) using BWA')
	parser_export.add_argument('restfile',
		type=str,
		help='Numpy compressed file containing restriction site coordinates')
	parser_export.add_argument('plotfile',
		type=str,
		help='Numpy compressed file containing restriction site coordinates')
	parser_export.set_defaults(func=exportToPlot)

	#
	parser_export = subparsers.add_parser('markdup',
		description='Add duplicate info for exported results')
	parser_export.add_argument('inifile',
		type=str,
		help=descIniFile)
	parser_export.add_argument('pdframe',
		type=str,
		help='The main file previously exported')
	parser_export.add_argument('extra',
		type=str,
		help='The _extra file created during previous export')
	parser_export.add_argument('outfile',
		type=str,
		help='New export file')
	parser_export.set_defaults(func=markDuplicates)


	args = parser.parse_args(sys.argv[1:])
	log.printArgs(args)
	args.func(args)

if __name__ == '__main__':

	main()
