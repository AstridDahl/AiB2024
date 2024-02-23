import numpy as np
from timeit import default_timer as timer
from time import sleep
import matplotlib.pyplot as plt

# How to track time of one alignment. 
then = timer()
align(seq1, seq2, substitution_matrix, 10, 5)
now = timer()

running_time = now - then
print(running_time)

ns = # store lengths of the longest sequence to align in every case. 
print(ns)
y = [] # To store experimental running times divided by theoretical running times. 

# Experimental running time. 
for n in ns:
    then = timer()
        # insert align() and track time of the different test cases. 
    now = timer()
running_times = now - then

# ys.append(running_times/)



# Theoretical running time. 
for n in ns: 
    then = timer()
        # use sleep to make theoretical running times. 
    now = timer()
theoretical_running_times = now - then