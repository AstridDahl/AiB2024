
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


def pick_max(eop_list):
    m_odd_even, indexes_oe = match_odd_to_even(eop_list)[0], match_odd_to_even(eop_list)[1]
    m_even_odd, indexes_eo = match_even_to_odd(eop_list)[0], match_even_to_odd(eop_list)[1]

    if m_even_odd > m_odd_even:
        return m_even_odd, "e2o", indexes_eo
    elif m_odd_even > m_even_odd:
        return m_odd_even, "o2e", indexes_oe
    elif m_odd_even == m_even_odd:
        return m_odd_even, "same matchscore", indexes_eo


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
        first_loop = 'l' + int((diff/2)-1) * "f" + "rr" + (int(diff/2)-1) * "f"

    return first_loop


def fold_between_matches_left_right_S1(eop): 

    # List of S1 matches
    indexes_dict = pick_max(eop)[2]
    S1_list = list(indexes_dict.keys())
    sorted_S1_list = sorted(S1_list)

    # S1 
    S1 = ""

    for i in range(len(sorted_S1_list)- 1):

        if len(S1) == 0 and sorted_S1_list[i] + 2 == sorted_S1_list[i+1]:
            S1 += 'ff'
            
        elif len(S1) == 0 and sorted_S1_list[i] + 2 < sorted_S1_list[i+1]:
            diff = sorted_S1_list[i+1] - sorted_S1_list[i] - 1
            f_number = (diff - 3) // 2
            fs = "f" * f_number
            S1 += "fl" + fs + "rr" + fs 

        elif sorted_S1_list[i] + 2 == sorted_S1_list[i+1]:
            
            if len(S1) != 0 and sorted_S1_list[i-1] < sorted_S1_list[i] - 2:
                S1 += "lf" 
            elif len(S1) != 0 and sorted_S1_list[i-1] == sorted_S1_list[i] - 2:
                S1 += "ff"
            
        elif sorted_S1_list[i] + 2 < sorted_S1_list[i+1]:

            diff = sorted_S1_list[i+1] - sorted_S1_list[i] - 1
            f_number = (diff - 3) // 2
            fs = "f" * f_number

            if S1[-1] == "f":
                S1 += "fl" + fs + "rr" + fs  

            elif S1[-1] == "r":
                S1 += "ll" + fs + "rr" + fs  
             
    return S1


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


def beginning_and_end(eop):
    indexes_dict = pick_max(eop)[2]
    min_key = min(indexes_dict.keys())
    min_value = indexes_dict[min_key]

    add_beginning = min_key - 0
    add_end = len(eop) - min_value - 1

    return add_beginning, add_end

def full_fold(eop):
    middle = first_fold(eop)
    left = fold_between_matches_left_right_S1(eop)
    right = fold_between_matches_right_left_S2(eop)
    start_and_end = beginning_and_end(eop)
    start = "f" * start_and_end[0]
    end = "f" * start_and_end[1]
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
S7 = even_odd("pphpphhpppphhpppphhpppphh")
print("full fold S7", full_fold(S7))
S8 = even_odd("ppphhpphhppppphhhhhhhpphhpppphhpphpp")
print("full fold S8", full_fold(S8))
S9 = even_odd("pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh")
print("full fold S9", full_fold(S9))
S10 = even_odd("hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh")
print("full fold S10", full_fold(S10))
S11 = even_odd("pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp")
print("full fold S11", full_fold(S11))
S12 = even_odd("hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh")
print("full fold S12", full_fold(S12))
S13 = even_odd("hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph")
print("full fold S13", full_fold(S13))
S14 = even_odd("pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh")
print("full fold S14", full_fold(S14))
S15 = even_odd("ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh")
print("full fold S15", full_fold(S15))

## --------------------- SCORES ------------------------- ##
'''
S1:
python hpview3k.py hhppppphhppphppphp fflfrrflfrrflrrlf
My score: 2
- opt: 4

S2:
python hpview3k.py hphphhhppphhhhpphh fffffffrrfffflrrl
My score: 5
- opt: 8

S3:
python hpview3k.py phpphphhhphhphhhhh ffffffffffrrfffff
My score: 4
- opt: 9

S4:
python hpview3k.py hphpphhphpphphhpphph ffflrrlffrrfflrrlff
My score: 4
- opt: 9

S5:
python hpview3k.py hhhpphphphpphphphpph fflrrlfffrrffffffff
My score: 4
-10

S6:
python hpview3k.py hhpphpphpphpphpphpphpphh flrrllfrrflrrlfrrfllrrl
My score : 6
-9

S7: 
python hpview3k.py pphpphhpppphhpppphhpppphh ffflrrlffffrrfffflfrrflf
My score: 2
-8

S8:
python hpview3k.py ppphhpphhppppphhhhhhhpphhpppphhpphpp ffffflrrllfrrflffrrlrrllfrrfllrrlff
My score: 4
-14

S9:
python hpview3k.py pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh ffflrrllrrllfrrflfffffrrflffrrffllrrllfrrflffff
My score: 11
-23

S10:
python hpview3k.py hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh flfffrrffflfffflrrllrrlfrrflrrllrrlfffflfffrrfffl
My score: 10
-21

S11: 
python hpview3k.py pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp ffffffffffffflfrrflfffffffffrrffffffffffffflfrrflffffffffff 
My score: 20
-36

S12:
python hpview3k.py hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh fffffffffffffffflrrllrrllfrrflfrrflfrrfllrrllrrlfffffffffffffff
My score: 18
-42

S13:
python hpview3k.py hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph fffflfrrflfffffffffflffrrfflfffffffffflrrlrrfffffffflrrlfffffffffflffrrffllrrllrrlff
My score: 22
-53

S14: 
python hpview3k.py pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh fffffffffflfrrflfffffffffflfrrfflrrllrrlfflrrlfffrrfflrrlfffffffflffffrrfffflfffffflfrrfllffrrfflfff
My score: Wrong length
-48

S15: 
python hpview3k.py ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh fffflrrlfflrrlfffflrrlfffflfffrrffflfffflrrlfffffffrrffffffflrrlfffffffflrrlfffflfrrflfffflffrrfflf
My score: 29
-50

'''