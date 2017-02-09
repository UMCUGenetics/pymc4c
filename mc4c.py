import argparse
import sys
import os

import log
import prep

import mc4ctools as mc


def makePrimerFasta(args):
	settings = mc.loadIni(args.inifile)
	primerSeqs = mc.getPrimerSeqs(settings)
	mc.writePrimerFasta(primerSeqs, args.outfile)


def cleaveReads(args):
	settings = mc.loadIni(args.inifile)
	primerLens = [len(x) for x in settings['prm_seq']]
	prmCuts = mc.combinePrimers(args.bamfile,primerLens)
	mc.applyCuts(args.fastqfile,args.outfile,prmCuts)


def splitReads(args):
	settings = mc.loadIni(args.inifile)
	restSeqs = settings['re_seq']
	# TODO: Substitute reference genome with reads (?)
	mc.findRestrictionSeqs(args.fastqfile,args.outfile,restSeqs)


# Huge wall of argparse text starts here
def main():
	descIniFile = 'File containing experiment specific details'
	descFqFile = 'Fastq file containing actual data from sequencing'

	parser = argparse.ArgumentParser(
		description="TODO: Replace this description")
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
		help='Bam file after makeprimerfa results were mapped by bowtie2') # then sorted and indexed using samtools'?
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


	args = parser.parse_args(sys.argv[1:])
	log.printArgs(args)
	args.func(args)

if __name__ == '__main__':
	main()
