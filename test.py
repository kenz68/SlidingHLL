import random

#paste this test below shll class.

def CountExact(Arr):
    cnt = 0
    n = len(Arr)
    for i in range(n):
        if Arr[i] == 1:
            cnt = cnt + 1
    return cnt


H = [0 for i in range(1000000)]
t = 0
dt = 0.01

A = SlidingHyperLogLog(10,10000)

n = random.randint(50000,60000)


for i in range(n):
    num = random.randint(0,999999)
    A.Add(num,t)
    H[num] = 1
    t = t + dt

ne = CountExact(H)
E = A.EstimateCardinality(t,t)
error = math.fabs((ne-E))/ne * 100
print("n = " + str(n))
print("Actual count: " + str(ne))
print("Estimated count: " + str(E))
print("Percentage of error: " + str(round(error,2)) + "%")






