import argparse
import sys
import os

import log
import prep

import rfs
import m2p

def main():
	parser = argparse.ArgumentParser(
		description="TODO: Replace this description")
	subparsers = parser.add_subparsers()

	parser_listIndex = subparsers.add_parser('listindex',
		description='Provide a list of all indexes in the matlab file')
	parser_listIndex.add_argument('infile',
		type=str,
		help='Matlab data object file')
	parser_listIndex.set_defaults(func=m2p.listIndexes)

	# File conversion
	parser_s01 = subparsers.add_parser('s01',
		description='S01 Preprocessing Data Information')
	parser_s01.add_argument('infile',
		type=str,
		help='Tab separated value file for conversion')
	parser_s01.add_argument('-soibeds',
		type=str, nargs='*',
		help='Bed files containing ...')
	parser_s01.add_argument('-antbeds',
		type=str, nargs='*',
		help='Bed files containing annotations ...')
	parser_s01.set_defaults(func=prep.prepareMeta)

	#
	parser_rfs01 = subparsers.add_parser('makeprimerfa',
		description='rfs Align Primers')
	parser_rfs01.add_argument('infile', # data_info
		type=str,
		help='Matlab data object file')
	parser_rfs01.add_argument('output', # rfs_prm
		type=str,
		help='Fasta file with primer sequences')
	parser_rfs01.add_argument('id',
		type=str, default=None,
		help='Id of the primer used, matching id column in infile')
	parser_rfs01.set_defaults(func=rfs.makePrimerFasta)

	#
	parser_rfs02 = subparsers.add_parser('cleavereads',
		description='Cleaving reads')
	parser_rfs02.add_argument('info_file', # data_info
		type=str,
		help='Matlab data object file')
	parser_rfs02.add_argument('rfs_file',
		type=str,
		help='Bam file after makeprimerfa results were mapped by bowtie2, then sorted and indexed using samtools')
	parser_rfs02.add_argument('infasta',
		type=str,
		help='Fasta file containing actual data from sequencing') # cmb_file
	parser_rfs02.add_argument('outfasta',
		type=str,
		help='Fasta file to dump cleaved sequences into')
	parser_rfs02.add_argument('id',
		type=str, default=None,
		help='Id of the primer used, matching id column in infile') # This id should be removed at some point later on...
	parser_rfs02.set_defaults(func=rfs.cleaveReads)

	#
	parser_gen = subparsers.add_parser('splitreads',
		description='Cleaving reads')
	parser_gen.add_argument('info_file', # data_info
		type=str,
		help='Matlab data object file')
	parser_gen.add_argument('infasta',
		type=str,
		help='Fasta file containing actual data from sequencing') # cmb_file
	parser_gen.add_argument('outfasta',
		type=str,
		help='Fasta file to dump split sequences into')
	parser_gen.add_argument('id',
		type=str, default=None,
		help='Id of the primer used, matching id column in infile') # This id should be removed at some point later on...
	parser_gen.set_defaults(func=rfs.splitReads)


	args = parser.parse_args(sys.argv[1:])
	log.printArgs(args)
	args.func(args)

if __name__ == '__main__':
	main()
