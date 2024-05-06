# ------ PACKAGES ------ #
from Bio import Phylo
import numpy as np
import re
import io
#------------------------#

########################
##       AIB P6       ##
########################

'''
This project is about protein folding in the 2D HP Model.
The objective is to fold a set of HP strings.

You can:

(1) Implement the 1/4 approximation algorithm for folding a hp-string in the 2D HP Model. 
You program should take a hp-string as input and output the energy of the computed fold for this string. 
It should also be possible to output a representation of the computed fold in a format that can be used 
by hpview3k.py (or hpview.py if you are using Python 2), see below, such that its energy can be validated.

(2) Implement another approach for finding good (hopefully near-optimal) folds in the 2D HP Model. 
The approach can for example be an algorithm extending or refining the folds found by 1/4 approximation
algorithm, or it can also be any other heuristic that you can think of or find in the litterature. 
All that is required is that you can explain your approach such that it is reproducible by others. 
Your approach should take a hp-string as input and output the energy of the best found fold for this string. 
It should also be possible to output a representation of the best found fold in a format that can be used 
by hpview.py (see below) such that its energy can be validated.

- hpview3k.py
Print the fold of a hp-string in ascii-format together with its score. 
TA folding in relative format is a sequence over f,l, and r: 
which describes  (f)orward, (l)left, or (r)ight. 
A folding in absolute format is a sequence over n,s,e, and w:
which describes (n)orth, (s)outh, (e)east, or (w)est. 

For example:
$ python hpview3k.py hphhhhhh ffffrrf 

h - p - h - H - h  
        *   *   |  
.   .   H - h - H  
                   
Score: 2

$ python hpview3k.py hphhhhhh eeeesww

h - p - h - H - h  
        *   *   |  
.   .   H - h - H  
                   
Score: 2

$ python hpview3k.py hphhhhhh eeeeswn

h - p - h - X - h  
            |   |  
.   .   .   h - H  
                   
Illegal fold after 7 steps

Examples:
Below are 15 hp-strings and the free energy of their optimal folds. 

 1: hhppppphhppphppphp -4
 2: hphphhhppphhhhpphh -8
 3: phpphphhhphhphhhhh -9
 4: hphpphhphpphphhpphph -9
 5: hhhpphphphpphphphpph -10
 6: hhpphpphpphpphpphpphpphh -9
 7: pphpphhpppphhpppphhpppphh -8
 8: ppphhpphhppppphhhhhhhpphhpppphhpphpp -14
 9: pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh -23
10: hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh -21
11: pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp -36
12: hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh -42
13: hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph -53
14: pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh -48
15: ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh -50

Presentation:
Your presentation must include (1)  a presentation of the folding method(s) that you have implemented, 
and (2) the best fold, and the time it took to find it, for each of the above 15 benchmark strings. 
'''

str_1 = "hhppppphhppphppphp"
str_2 = "hphphhhppphhhhpphh"

def even_odd(str):
    list_hp = list(str)
    list_even_odd_p = []
    for i in range(len(list_hp)):
        if i % 2 == 1 and list_hp[i] == 'h': # Odd and 'h'
            list_even_odd_p.append('o')
        elif i % 2 == 0 and list_hp[i] == 'h':
            list_even_odd_p.append('e')
        else:
            list_even_odd_p.append('p')
    return list_even_odd_p

print("str_1", even_odd(str_1))
print("str_2", even_odd(str_2))

eop_1 = ['e', 'o', 'p', 'p', 'p', 'p', 'p', 'o', 'e', 'p', 'p', 'p', 'e', 'p', 'p', 'p', 'e', 'p']
eop_2 = ['e', 'p', 'e', 'p', 'e', 'o', 'e', 'p', 'p', 'p', 'e', 'o', 'e', 'o', 'p', 'p', 'e', 'o']

def match_even_to_odd(eop_list):
    match = 0
    i = 0
    j = len(eop_list) - 1
    limit = int(1/2*len(eop_list))
    indexes_matches = {}
    while i <= limit and limit <=j:
        if eop_list[i] != 'e':
            i = i + 1
        elif eop_list[i] == 'e':
            if eop_list[j] != 'o':
                j = j - 1
            elif eop_list[j] == 'o':
                indexes_matches[i] = j
                match = match +1 
                i = i + 1
                j = j - 1
    return match, indexes_matches
     
print("matches eop_1 even_odd", match_even_to_odd(eop_1))
print("matches eop_2 even_odd", match_even_to_odd(eop_2))

def match_odd_to_even(eop_list):
    match = 0
    i = 0
    j = len(eop_list) - 1
    limit = int(1/2*len(eop_list))
    indexes_matches = {}
    while i <= limit and limit <=j:
        if eop_list[i] != 'o':
            i = i + 1
        elif eop_list[i] == 'o':
            if eop_list[j] != 'e':
                j = j - 1
            elif eop_list[j] == 'e':
                indexes_matches[i] = j
                match = match +1
                i = i + 1
                j = j - 1
    return match, indexes_matches

print("matches eop_1 odd_even", match_odd_to_even(eop_1))
print("matches eop_2 odd_even", match_odd_to_even(eop_2))

def pick_max(eop_list):
    m_odd_even, indexes_oe = match_odd_to_even(eop_list)[0], match_odd_to_even(eop_list)[1]
    m_even_odd, indexes_eo = match_even_to_odd(eop_list)[0], match_even_to_odd(eop_list)[1]

    if m_even_odd > m_odd_even:
        return m_even_odd, "e2o", indexes_eo
    elif m_odd_even > m_even_odd:
        return m_odd_even, "o2e", indexes_oe
    elif m_odd_even == m_even_odd:
        return m_odd_even, "same matchscore"


print("max matches eop_1", pick_max(eop_1))
print("max matches eop_2", pick_max(eop_2))

def first_fold(eop):
    indexes_dict = pick_max(eop)[2]
    
    # Get largest key in dict
    # The last match-pair before first fold
    max_key = max(indexes_dict.keys())
    max_value = indexes_dict[max_key]

    # How many aa between them
    diff = max_value - max_key - 1

    S1_list = list(indexes_dict.keys())
    sorted_S1_list = sorted(S1_list, reverse= True)
    first_element = sorted_S1_list[0]
    second_element = sorted_S1_list[1]

    if 2 + second_element == first_element:
        first_loop = int(diff/2) * "f" + "rr" + (int(diff/2)-1) * "f"
    else:
        first_loop = 'l' + int((diff/2)/2) * "f" + "rr" + (int(diff/2)-1) * "f"

    return first_loop

print("first fold eop_1", first_fold(eop_1))
print("first fold eop_2", first_fold(eop_2))


def fold_between_matches_left_right_S1(eop): # S1

    # List of S1 matches
    indexes_dict = pick_max(eop)[2]
    S1_list = list(indexes_dict.keys())
    sorted_S1_list = sorted(S1_list)

    # S1 
    S1 = ""

    for i in range(len(sorted_S1_list)- 1):

        if sorted_S1_list[i] + 2 == sorted_S1_list[i+1]:
            S1 += "ff"

        elif sorted_S1_list[i] + 2 < sorted_S1_list[i+1]:
            diff = sorted_S1_list[i+1] - sorted_S1_list[i] - 1
            f_number = (diff - 3) // 2
            fs = "f" * f_number
            S1 += "fl" + fs + "rr" + fs     
    return S1

print("between S1 matches eop_1", fold_between_matches_left_right_S1(eop_1))
print("between S1 matches eop_2", fold_between_matches_left_right_S1(eop_2))

def fold_between_matches_right_left_S2(eop):

    # List of S2 matches
    indexes_dict = pick_max(eop)[2]
    values_list = list(indexes_dict.values())
    sorted_S2_list = sorted(values_list)  

    # S2
    S2 = ""

    for i in range(len(sorted_S2_list)- 1):

        if sorted_S2_list[i] + 2 == sorted_S2_list[i+1]:
            S2 += "ff"

        elif sorted_S2_list[i] + 2 < sorted_S2_list[i+1]:
            diff = sorted_S2_list[i+1] - sorted_S2_list[i] - 1
            f_number = (diff - 3) // 2
            fs = "f" * f_number
            S2 += "l" + fs + "rr" + fs + "l" 

    return S2

print("between S2 matches eop_1", fold_between_matches_right_left_S2(eop_1))
print("between S2 matches eop_2", fold_between_matches_right_left_S2(eop_2))

def beginning_and_end(eop):
    indexes_dict = pick_max(eop)[2]
    min_key = min(indexes_dict.keys())
    min_value = indexes_dict[min_key]
    diff = min_value - min_key
    len_eop = len(eop)
    missing = len_eop - diff - 1
    return missing

print("start and end", beginning_and_end(eop_1))
print("start and end", beginning_and_end(eop_2))

def full_fold(eop):
    middle = first_fold(eop)
    left = fold_between_matches_left_right_S1(eop)
    right = fold_between_matches_right_left_S2(eop)
    start_and_end = beginning_and_end(eop)
    start = "f" * int(start_and_end/2)
    end = "f" * int(start_and_end/2)
    return start+left+middle+right+end


## --------- TEST --------- ##
S1 = even_odd("hhppppphhppphppphp")
print("full fold S1", full_fold(S1))
S2 = even_odd("hphphhhppphhhhpphh")
print("full fold S2", full_fold(S2))
S3 = even_odd("phpphphhhphhphhhhh")
print("full fold S3", full_fold(S3))
S4 = even_odd("hphpphhphpphphhpphph")
print("full fold S4", full_fold(S4))
S5 = even_odd("hhhpphphphpphphphpph")
print("full fold S5", full_fold(S5))
S6 = even_odd("hhpphpphpphpphpphpphpphh")
print("full fold S6", full_fold(S6))

'''
S1:
python hpview3k.py hhppppphhppphppphp fflfrrflfrrflrrlf
My score: 2

S2:
python hpview3k.py hphphhhppphhhhpphh fffffffrrfffflrrl
My score: 5

S3:
python hpview3k.py phpphphhhphhphhhhh ffffffffrrfffffff
My score: 4

S4:
python hpview3k.py hphpphhphpphphhpphph ffflrrfffrrfflrrlff
- Illegal fold
-9

S5:
python hpview3k.py hhhpphphphpphphphpph ffflrrffffrrfffffff
-10

S6:
python hpview3k.py hhpphpphpphpphpphpphpphh flrrflfrrflrrlfrrfllrrl
-9

S7:
python hpview3k.py pphpphhpppphhpppphhpppphh -8

S8:
python hpview3k.py ppphhpphhppppphhhhhhhpphhpppphhpphpp -14

S9:
python hpview3k.py pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh -23

S10:
python hpview3k.py hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh -21

S11: 
python hpview3k.py pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp -36

S12:
python hpview3k.py hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh -42

S13:
python hpview3k.py hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph -53

S14: 
python hpview3k.py pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh -48

S15: 
python hpview3k.py ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh -50

'''