#!/usr/bin/env python
#
#  """HTSeqQA_mapper.py"""
  #
__author__  ="Alexey Siretskiy"
__email__   ="alexey.siretskiy@it.uu.se"
__date__    ="2014"
__doc__			="""This script is a part of the HTSeq-Hadoop package for
a  high-throughput sequencing reads or ``SAM`` files  analysis. The script takes
a file with high-throughput sequencing reads  from ``STDIN``, and emerges to ``STDOUT``
the bases and  base-call quality scores by position  within the reads.
Should be used with :mod:`HTSeq_Hadoop.HTSeqQA_reducer` script.
The input data at the moment is ``SAM``."""
#%s, IT dept, Uppsala University, Sweden, %s, %s"""%(__author__ ,__email__, __date__)


import sys, time, os.path, argparse
from itertools import *
import numpy
import HTSeq

def get_reads_length(reads, n_lines):
	""" calculate the read length from the data 
based on the first 100 lines of the input file
				
:param reads: file chunk to determine the read lenght.
:type name: `str`.
:param n_lines: the number of lines to analyze in the chunk
:type name: `int`
:return  readlen:  -- the number of lines.
:rtype: `int`
"""

	readlen = 0
	for r in islice( reads, n_lines ):
		if len( r ) > readlen:
			readlen = len( r )
	return readlen


def main():
	
	parser = argparse.ArgumentParser(	description=__doc__ )
	
	parser.add_argument( "-t", "--type", dest="type",\
			choices=["sam",  "fastq"],\
			default = "sam", help="type of read_file (one of: sam [default], fastq" )
	#parser.add_argument( "-r", "--readlength", type=int, dest="readlen",\
	#		help="the maximum read length (when not specified, the script guesses from the file" )
	parser.add_argument( "-n", "--nosplit", action="store_true", dest="nosplit",\
			help="do not split reads in unaligned and aligned ones" )
	parser.add_argument( "-m", "--maxqual", type=int, dest="maxqual", default=70,\
			help="the maximum quality score that appears in the data (default: 70)" )

	opts = parser.parse_args()


	if len( sys.argv ) < 2:
		print parser.parse_args(['-h'])

	"""read the STDIN"""
	readfilename = sys.stdin #args[0]

	if opts.type == "sam":
		readfile = HTSeq.SAM_Reader( readfilename )
#		print readfile
		isAlnmntFile = True
#	elif opts.type == "bam":
#		readfile = HTSeq.BAM_Reader( readfilename )
#		isAlnmntFile = True
#	elif opts.type == "solexa-export":
#		readfile = HTSeq.SolexaExportReader( readfilename )
#		isAlnmntFile = True
	elif opts.type == "fastq":
		readfile = HTSeq.FastqReader( readfilename )
		isAlnmntFile = False
#	elif opts.type == "solexa-fastq":
#		readfile = HTSeq.FastqReader( readfilename, "solexa" )
#		isAlnmntFile = False
	else:
		sys.error( "Oops." )

	twoColumns = isAlnmntFile and not opts.nosplit
	#print "twocolumns", twoColumns, 'file type', opts.type

# **** Get read length ****
	if isAlnmntFile:
		reads = ( a.read for a in readfile )
	else:
		reads = readfile
	
	""" use first 100 lines to obtain the read length"""	
	readlen = get_reads_length(reads, 100)	#print readlen
	max_qual = opts.maxqual 
	""" SAM files can have qualities higher than FASTQ files..."""


# **** Initialize count arrays ****

	base_arr_U = numpy.zeros( ( readlen, 5 ), numpy.int )
	qual_arr_U = numpy.zeros( ( readlen, max_qual+1 ), numpy.int )
	if twoColumns:
		base_arr_A = numpy.zeros( ( readlen, 5 ), numpy.int )
		qual_arr_A = numpy.zeros( ( readlen, max_qual+1 ), numpy.int )


# **** Main counting loop ****

	i = 0
	try:
		for a in readfile:
			if isAlnmntFile:
				r = a.read
			else:
				r = a
			if twoColumns and (isAlnmntFile and a.aligned):
				r.add_bases_to_count_array( base_arr_A )
				r.add_qual_to_count_array( qual_arr_A )
			else:
				r.add_bases_to_count_array( base_arr_U )
				r.add_qual_to_count_array( qual_arr_U )   
			i += 1
	except:
		sys.stderr.write( "Error occured in: %s\n" %
		readfile.get_line_number_string() )
		raise
	#print i, "reads processed"
	#print "number of aligned reads", base_arr_A[50,:].sum()
	#print "number of unaligned reads", base_arr_U[60,:].sum()
	
	#print 'length quals', len(qual_arr_U[0,:])
	#print 'lenght read', len(base_arr_A), base_arr_A.sum()
	
	
	for i  in range(len(base_arr_A)):
		print "%s\t%s" %( "base_arr_U_"+str(i), ' '.join(str(x) for x in base_arr_U[i]) )
		print "%s\t%s" %( "qual_arr_U_"+str(i), ' '.join(str(x) for x in qual_arr_U[i]) )
	print "%s\t%s" %("nreads_U_U", str(base_arr_U[0,:].sum()) )

	
	if twoColumns: 
		for i  in range(len(base_arr_A)):
			print "%s\t%s" %( "base_arr_A_"+str(i), ' '.join(str(x) for x in base_arr_A[i]) )
			print "%s\t%s" %( "qual_arr_A_"+str(i), ' '.join(str(x) for x in qual_arr_A[i]) )
		print "%s\t%s" %("nreads_A_A", str(base_arr_A[0,:].sum()) )
		

if __name__ == "__main__":
	main()
