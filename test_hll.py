# prompt: write test with random data that fake streaming data with time t increase as timestamp

import random
from hll import HyperLogLog
import sys
b = 10
n = 1000000
# Create a HyperLogLog object
hll = HyperLogLog(b)
actual_set = set()

# Generate random data and add it to the hll
for i in range(n):
    num = random.randint(0, n)
    actual_set.add(num)
    hll.Add(num)

est = hll.EstimateCardinality()

# Print the results
print("n =", len(actual_set))
print("size of actual set: ", sys.getsizeof(actual_set), " bytes")
print("Estimated cardinality:", est)
print("size of hll: ", sys.getsizeof(hll), " bytes")
