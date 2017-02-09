import argparse
import sys
import os

import log
import prep

import rfs
import m2p

def main():
	descIniFile = 'File containing experiment specific details'
	descFqFile = 'Fastq file containing actual data from sequencing'

	parser = argparse.ArgumentParser(
		description="TODO: Replace this description")
	subparsers = parser.add_subparsers()

	#
	parser_rfs01 = subparsers.add_parser('makeprimerfa',
		description='rfs Align Primers')
	parser_rfs01.add_argument('inifile', # data_info
		type=str,
		help=descIniFile)
	parser_rfs01.add_argument('outfile', # rfs_prm
		type=str,
		help='Fasta file with primer sequences')
	parser_rfs01.set_defaults(func=rfs.makePrimerFasta)

	#
	parser_rfs02 = subparsers.add_parser('cleavereads',
		description='Cleave reads by primer sequences')
	parser_rfs02.add_argument('inifile', # data_info
		type=str,
		help=descIniFile)
	parser_rfs02.add_argument('bamfile',
		type=str,
		help='Bam file after makeprimerfa results were mapped by bowtie2') # then sorted and indexed using samtools'?
	parser_rfs02.add_argument('fastqfile',
		type=str,
		help=descFqFile) # cmb_file
	parser_rfs02.add_argument('outfile',
		type=str,
		help='Fastq file to dump primer cleaved sequences into')
	parser_rfs02.set_defaults(func=rfs.cleaveReads)

	#
	parser_gen = subparsers.add_parser('splitreads',
		description='Split reads by restriction site sequences')
	parser_gen.add_argument('inifile', # data_info
		type=str,
		help=descIniFile)
	parser_gen.add_argument('fastqfile',
		type=str,
		help=descFqFile) # cmb_file
	parser_gen.add_argument('outfile',
		type=str,
		help='Fasta file to dump restriction site split sequences into')
	parser_gen.set_defaults(func=rfs.splitReads)


	args = parser.parse_args(sys.argv[1:])
	log.printArgs(args)
	args.func(args)

if __name__ == '__main__':
	main()
