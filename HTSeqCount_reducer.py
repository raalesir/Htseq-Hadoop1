#!/usr/bin/env python
"""A reducer for summing up the gene -- count lists collected from the mappers.
The input is read from ``STDIN``, the output is written to ``STDOUT``"""


from itertools import groupby
from operator import itemgetter
import sys
import argparse

#parser = argparse.ArgumentParser(description='Reducer for the  HTSeq-count utility')
#parser.add_argument('-o','--output', help='file to write the result to', required=True)

def read_mapper_output(file, separator='\t'):
	for line in file:
		yield line.rstrip().split(separator, 1)

def main(separator='\t'):
	data = read_mapper_output(sys.stdin, separator=separator)
	# groupby groups multiple word-count pairs by word,
	# and creates an iterator that returns consecutive keys and their group:
	#   current_word - string containing a word (the key)
	#   group - iterator yielding all ["&lt;current_word&gt;", "&lt;count&gt;"] items

	for current_word, group in groupby(data, itemgetter(0)):
		try:
			total_count = sum(int(count) for current_word, count in group)
			print "%s%s%d" % (current_word, separator, total_count)
		except ValueError:
		# count was not a number, so silently discard this item
			pass

if __name__ == "__main__":
	main()
