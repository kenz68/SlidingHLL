# prompt: write test with random data that fake streaming data with time t increase as timestamp

import random
import math
import time

from shll import SlidingHyperLogLog

# Define the parameters of the Sliding HyperLogLog
b = 3
W = 1

# Create a Sliding HyperLogLog object
shll = SlidingHyperLogLog(b, W)

# Generate random data and add it to the shll
t = 0
dt = 0.01
n = 8
for i in range(8):
    num = random.randint(0, 8)
    shll.Add(num, t)
    t += 1

print(shll.calculate_cardinality_buckets)
# Estimate the cardinality of the data
E = shll.EstimateCardinality(t, t)

# Print the results
print("n =", n)
print("Estimated cardinality:", E)
