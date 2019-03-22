"""RegExpManager module offers
convenience functions to check strings
against regular expression, e.g. to match...
- chromosomes
- file types
- file names
"""

import re as re

# TODO Evaluate raison d'etre
# So far, there is no reason why this should be a class/object

class RegExpManager(object):
	def __init__(self):
		# fix definitions for regexps are stored within
		# the functions since a RegExpManager is only needed
		# every now and then
		pass
	def get_chrom_matcher(self,c):
		"""get_chrom_matcher
		Receives a string that designates either a pre-defined
		set of chromosomes (like autosomes) or is interpreted
		as a compilable regexp string. Returns a function that
		expects another string to test (True/False testing)
		"""
		# pre-defined expressions
		chromdict = { "autosomal":r"\bchr[0-9]{1,2}[ab]?\b", # match single autosome like chr22, chr1, chr2a (chimp)
					"gonosomal":r"\bchr[XYZW]\b", # match single gonosome like chrX, chrW (e.g. chicken)
					"mitochondrial":r"\bchrM\b", # match single mitochondrial chromosome chrM
					"all":r"\bchr\w+\b" # match all types of chromosomes starting with chr
					}
		comp = None
		try:
			found = chromdict[c]
			comp = re.compile(found)
		except KeyError:
			comp = re.compile(c)
		finally:
			assert comp is not None, "get_chrom_matcher: could not compile expression {}".format(c)
			def matcher(to_test):
				mobj = comp.match(to_test)
				return mobj is not None
			return matcher
