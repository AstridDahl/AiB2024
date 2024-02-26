import numpy as np
from timeit import default_timer as timer
from time import sleep
import matplotlib.pyplot as plt

# from global_affine import align

# Compare running time with theoretical bounds
# Compare running time for alignment with linear gap cost to that of affine gap cost

# How to track time of one alignment. 
# then = timer()
# align(seq1, seq2, substitution_matrix, 10, 5)
# now = timer()

# running_time = now - then
# print(running_time)

# ns = # store lengths of the longest sequence to align in every case.
# print(ns)
# y = [] # To store experimental running times divided by theoretical running times. 

# # Experimental running time. 
# for n in ns:
#     then = timer()
#         # insert align() and track time of the different test cases. 
#     now = timer()
# running_times = now - then

# # ys.append(running_times/)

# # Theoretical running time. 
# for n in ns: 
#     then = timer()
#         # use sleep to make theoretical running times. 
#     now = timer()
# theoretical_running_times = now - then


# linear
# time = user + sys
# test1, 12: 0.142s

# test4, 200: 0.202s

# test5, 647: 0.570s

# test6, 1781: 3.423s


# affine
# time = user + sys
# test1, 12: 0.144s

# test4, 200: 0.416s

# test5, 647: 0.528s

# test6, 1781: 3.166s

ns = [12, 200, 647, 1781]
lin_run = [0.142, 0.202, 0.57, 3.423]
aff_run = [0.144, 0.416, 0.528, 3.166]
lin_ys = [0.142/12**2, 0.202/200**2, 0.57/647**2, 3.423/1781**2]
aff_ys = [0.144/12**3, 0.416/200**3, 0.528/647**3, 3.166/1781**3]


plt.plot(ns,lin_run)
plt.xlabel("Input size")
plt.ylabel("T(n)")
plt.savefig("linT_plot.pdf")

plt.plot(ns,aff_run)
plt.xlabel("Input size")
plt.ylabel("T(n)")
plt.savefig("affT_plot.pdf")

plt.plot(ns,lin_ys)
plt.xlabel("Input size")
plt.ylabel("T(n) / n^2")
plt.savefig("lin_plot.pdf")

plt.plot(ns,lin_ys)
plt.xlabel("Input size")
plt.ylabel("T(n) / n^3")
plt.savefig("aff_plot.pdf")
