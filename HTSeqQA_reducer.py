#!/usr/bin/env python
#
#
#  """HTSeqQA_reducer.py"""
#
#
__author__  ="Alexey Siretskiy"
__email__   ="alexey.siretskiy@it.uu.se"
__date__    ="2014"
__doc__ 		="""This script is a part of the HTSeq-Hadoo package for analyzing   high-throughput sequencing reads or ``SAM`` files. The script  reads the data from ``STDIN`` and producing plots showing the distribution of called bases and  base-call quality scores by position  within the reads. The plots are sent to ``STDOUT`` in the ``SVG`` format. Should be used with :mod:`HTSeq_Hadoop.HTSeqQA_mapper` script."""
#%s, IT dept, Uppsala Univervity, Sweden, %s, %s""" %(__author__ ,__email__, __date__)

import sys, time, os.path, argparse
from itertools import *
import numpy as np
import HTSeq
from operator import itemgetter

def read_mapper_output(file, separator='\t'):
	"""
			reads the STDIN lazily and ouputs the **<key,value>** separates by TAB
	"""
	for line in file:
		yield line.rstrip().split(separator, 1)



try:
	import matplotlib
except ImportError:
	sys.stderr.write("This script needs the 'matplotlib' library, which ")
	sys.stderr.write("was not found. Please install it." )
matplotlib.use('svg')
from matplotlib import pyplot



# **** Normalize result ****

def norm_by_pos( arr ):
	"""
	normalizes array by dividing each element in the row by the sum of the row 
	
	:param arr:  2D array 
	:type arr: `numpy.int`
	:return arr:   normalized array
	:rtype: `numpy.float`
	"""
	arr = np.array( arr, np.float )
	arr_n = ( arr.T / arr.sum( 1 ) ).T
	arr_n[ arr == 0 ] = 0
	return arr_n


def norm_by_start( arr ):
	"""
	normalizes array by the first element of the array, resulted after summing up elements in the rows
	
	:param arr:  2D array 
	:type arr: `numpy.int`
	:return arr: arr --  normalized array
	:rtype: `numpy.float`
	"""
	arr = np.array( arr, np.float )
	arr_n = ( arr.T / arr.sum( 1 )[ 0 ] ).T
	arr_n[ arr == 0 ] = 0
	return arr_n


def plot_bases( arr ):
	"""
		prepares the Matplotlib canvas
	"""
	readlen = len(arr)
	print readlen
	xg = np.arange( readlen )   
	pyplot.plot( xg, arr[ : , 0 ], marker='.', color='red')
	pyplot.plot( xg, arr[ : , 1 ], marker='.', color='darkgreen')
	pyplot.plot( xg, arr[ : , 2 ], marker='.',color='lightgreen')
	pyplot.plot( xg, arr[ : , 3 ], marker='.',color='orange')
	pyplot.plot( xg, arr[ : , 4 ], marker='.',color='grey')
	pyplot.axis( (0, readlen-1, 0, 1 ) )
	pyplot.text( readlen*.70, .9, "A", color="red" )
	pyplot.text( readlen*.75, .9, "C", color="darkgreen" )
	pyplot.text( readlen*.80, .9, "G", color="lightgreen" )
	pyplot.text( readlen*.85, .9, "T", color="orange" )
	pyplot.text( readlen*.90, .9, "N", color="grey" )

	
def main(separator='\t'):

	# **** Parse command line ****

	parser = argparse.ArgumentParser(description=__doc__)	
	parser.add_argument( "-g", "--gamma", type=float, dest="gamma", default = 0.3,\
			help="the gamma factor for the contrast adjustment of the quality score plot, default 0.3" )
	opts = parser.parse_args()

	gamma = opts.gamma
	twoColumns = True
	
# input comes from STDIN (standard input)
	data = read_mapper_output(sys.stdin, separator=separator)
	# groupby groups multiple word-count pairs by word,
	# and creates an iterator that returns consecutive keys and their group:
	# current_word - string containing a word (the key)
	# group - iterator yielding all ["&lt;current_word&gt;","&lt;count&gt;"] items
	former_element = None
	array_types = {}
	for current_element, group in groupby(data, itemgetter(0)):
		try:
			total_count = sum(np.fromstring(line, sep=' ', dtype=int) for current_element, line in group)
		#	print "%s%s%s" % (current_element, separator, ' '.join(str(x) for x in total_count))
		#	print current_element, total_count
			if not former_element: 
				arr = total_count
				former_element = current_element
			else:
				if current_element.split('_')[:-1] == former_element.split('_')[:-1]:
					arr = np.vstack([arr,total_count])
				else:
				#print 'change array', former_element, current_element
					array_types['_'.join(former_element.split('_')[:-1])] = arr		
					former_element = current_element
					arr = total_count
		except ValueError:
		# count was not a number, so silently discard this item
			pass
	array_types['_'.join(former_element.split('_')[:-1])] = arr
	

	base_arr_U_n = norm_by_pos( array_types['base_arr_U'] )
	qual_arr_U_n = norm_by_start( array_types['qual_arr_U'] )
	nreads_U = array_types['nreads_U'] # base_arr_U[0,:].sum()
	if twoColumns:
		base_arr_A_n = norm_by_pos( array_types['base_arr_A'] )
		qual_arr_A_n = norm_by_start( array_types['qual_arr_A'] )
		nreads_A = array_types['nreads_A'] #base_arr_A[0,:].sum()
	
	# get the read length from the data
	readlen = len(base_arr_A_n)
	#get the max qual from the data
	max_qual = len(qual_arr_U_n[0,:])
#	print readlen, max_qual, len(base_arr_U_n)
#	sys.exit()
#	print "number of aligned reads", base_arr_A[50,:].sum()
#	print "number of unaligned reads", base_arr_U[60,:].sum()
#
#
## **** Make plot ****
#

	pyplot.figure()
	pyplot.subplots_adjust( top=.85 )
	pyplot.suptitle( os.path.basename("readfilename"), fontweight='bold' )

	if twoColumns:

		pyplot.subplot( 221 )
		plot_bases( base_arr_U_n )
		pyplot.ylabel( "proportion of base" )
		pyplot.title( "non-aligned reads\n%.0f%% (%.3f million)" %( 100. * nreads_U / (nreads_U+nreads_A), nreads_U / 1e6 ) )

		pyplot.subplot( 222 )
		plot_bases( base_arr_A_n )
		pyplot.title( "aligned reads\n%.0f%% (%.3f million)" %( 100. * nreads_A / (nreads_U+nreads_A), nreads_A / 1e6 ) )

		pyplot.subplot( 223 )
		pyplot.pcolor( qual_arr_U_n.T ** gamma, cmap=pyplot.cm.Greens,
		norm=pyplot.normalize( 0, 1 ) )
		pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
		pyplot.xlabel( "position in read" )
		pyplot.ylabel( "base-call quality score" )

		pyplot.subplot( 224 )
		pyplot.pcolor( qual_arr_A_n.T ** gamma, cmap=pyplot.cm.Greens,
		norm=pyplot.normalize( 0, 1 ) )
		pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
		pyplot.xlabel( "position in read" )

	else:

		pyplot.subplot( 211 )
		plot_bases( base_arr_U_n )
		pyplot.ylabel( "proportion of base" )
		pyplot.title( "%.3f million reads" % ( nreads_U / 1e6 ) )

		pyplot.subplot( 212 )
		pyplot.pcolor( qual_arr_U_n.T ** gamma, cmap=pyplot.cm.Greens,
		norm=pyplot.normalize( 0, 1 ) )
		pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
		pyplot.xlabel( "position in read" )
		pyplot.ylabel( "base-call quality score" )

	pyplot.savefig( sys.stdout )

if __name__ == "__main__":
	main()
