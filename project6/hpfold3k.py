#  hpfold3k.py
#
#  Implements the simple 1/4-approximation algorithm by Hart and Istrail
#
#  Usage:
#
#  hpfold.py <seq>
#
#  where <seq> is a string over the alphabet {h,p}. The output is a fold
#  in relative format.
#
#  A folding in relative format is a sequence over {f,l,r}, which
#  describes the fold as a sequence of steps (f)orward, (l)left, or
#  (r)ight.
#
#  History:
# 
#  26-Nov-2009: Initial version.
#  05-Dec-2015: Small changes (print and has_key) to adapt to Python 3.
#
#  Christian Storm Pedersen <cstorm@birc.au.dk>

import sys
from hpview3k import HPFold

#####################################################################
# Main
#####################################################################

def make_loop(k):
	if k == 2:
		loop = 'ff'
	else:
		loop = 'l' + ((k - 4) // 2) * 'f' + 'rr' + ((k - 4) // 2) * 'f' + 'l' # a gap of four means only the short terms should be added, 'lrrl'
	return loop
		
if __name__ == "__main__":
	
	seq = sys.argv[1]

	# Make lists of odd and even h's
	odd = [i for i in range(1, len(seq), 2) if seq[i] == 'h']
	even = [i for i in range(0, len(seq), 2) if seq[i] == 'h']
	
	# Make pairings of odd and even h's
	odd_first_pair = [(i,j) for (i,j) in zip(odd, reversed(even)) if j-i > 1]
	even_first_pair = [(i,j) for (i,j) in zip(even, reversed(odd)) if j-i > 1]

	print(odd_first_pair)
	print(even_first_pair)

	# Select largest pairing
	if len(odd_first_pair) > len(even_first_pair):
		pair = odd_first_pair
	else:
		pair = even_first_pair
	
	if len(pair) == 0:
		# No pairs of h's to match. Fold as a straight line
		fold = (len(seq) - 1) * 'f'
	else:
		# Make fold according to 1/4-approximation algorithm

		# Fold upper part
		fold = pair[0][0] * 'f'   		# how many forward we have to go to get to the first match of hs
		for i in range(1, len(pair)):
			fold = fold + make_loop(pair[i][0] - pair[i-1][0]) # the loops between all pairs
	
		# Make bend
		k = (pair[-1][1] - pair[-1][0]) // 2       # last pair is at index -1. k is how many forward we have to go before the bend, 'rr' is the actual bend
		fold = fold + k * 'f' + 'rr' + (k-1) * 'f' # (k-1) bc we already took the first step by the second r in 'rr'.
	
		# Fold lower part
		for i in range(1, len(pair)):
			fold = fold + make_loop(pair[-(i+1)][1] - pair[-i][1]) # same as before but we start from the end pair
		fold = fold + (len(seq) - pair[0][1] - 1) * 'f' # the final forwards to get to the start again
		
	# Print fold using HPFold class from hpview.py
	s = HPFold(seq)
	s.SetRelFold(fold)
	s.PrintFold()
