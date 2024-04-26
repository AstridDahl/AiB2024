# Measure running time
from timeit import default_timer as timer
from time import sleep

# Iterate over files in directory
import glob

# Running time for 14 distance matrices. 
d = dict()
for filepath in glob.iglob('unique_distance_matrices\*.phy'):
    then = timer()
    print(NJ_saitou_and_nei(filepath)) # Insert algorithm. 
    now = timer()
    running_time = now - then # Running time in seconds.
    print(running_time)
    d[filepath] = running_time

print(d) # Dictionary with filepath as key and running time as value. 

