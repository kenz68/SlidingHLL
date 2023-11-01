
import math
import heapq
#can import sha256 instead of sha1 too, but it's a bit slower!
from hashlib import sha1


def calculate_alpha_m(b):
    # the constant alpha_m in the article. m = 2^b, here we send b itself.
    # we calculate it strictly without multiplying with anything else.
    if (b > 16 or b < 4):
        raise ValueError("b=%d is not in range [4,16]" % b)
    if b == 4:
        return 0.673
    if b == 5:
        return 0.697
    if b == 6:
        return 0.709
    # for values b >=7, this equation is defined in the article.
    # 1 << b is the bit 1 that is shifted b bits to the left. so 1 << b = 2^b = m
    return 0.7213 / (1.0 + 1.079 / (1 << b))
    
    
def calculate_p_w(w, max_len):
    """
    Defined in the article:
    w = x(b+1),x(b+2),...
    p(w) = the position of the leftmost 1-bit of word w when written in binary form.
    we count from the left. e.g: p(0001) = 4, p (1...) = 1, p(0^k) = k+1
    NOTE: max_len is NOT the length of the type of w! (64/32 bits), because we truncated the first b bits for hashing!
    max_len is in fact the word length (memory representation) - b bits truncated
    """
    p_w = max_len - (w.bit_length()) + 1
    # p_w should be strictly in range 1...k+1
    if p_w <= 0:
        raise ValueError("w overflow")
    return p_w
    


class SlidingHyperLogLog(object):
    """ Implementation of Sliding HyperLogLog class.
    Properties of the class:
    m = 2^b - the number of registers
    b -  (log 2 of m above)
    alpha_m - the const used for corrction of hash bias (read article)
    W - maximum time window size
    LFPM - list of future possible maxima: a list of pairs (ti,p(wi))
    """


    def __init__(self, param, W):
        """
        The constructor of the class.
        We enable the end-user to use the constructor in one of two ways:
        1. Sending the log of the number of required registers directly: 
        param needs to be of type int, b = log_2(m) with b in range [4..16]
        2. Sending the allowed Std. error as a fraction (0 < r < 1)
        param then is of type float!
        Std. error is defined as in the article: (E' - E) / E where E is actual cardinality, E' is our estimate
        """

        if not type(W) == int:
            raise TypeError("Max. window size should be an integer")
        if W <= 0 :
            raise TypeError("Max. window size should be positive")
        self.W = W
        
        if type(param) == int:
        # here param is m, of type int.
        # here we choose not to send value error, rather to auto-correct inappropriate b values later
            b = param
        
        else:
            # we assume that param is std. error, of type float
            if not (0 < param < 1):
                raise ValueError("Std. error must be between 0 and 1.")
            # std. error = 1.04 / sqrt(m)
            b = int(math.ceil(math.log((1.04 / param) ** 2, 2)))
        
        # make adjustments for value of b, should not exceed range 4...16 (for both cases)
        if (b > 16):
            b = 16
        if (b < 4):
            b = 4

        #common code for both cases: set b, m, alpha_m, LFPM
        self.b = b
        # m = 2 ** b (2^b)
        self.m = 1 << b 
        self.alpha_m = calculate_alpha_m(b)
        # init. an empty list. LFPM is a list of pairs
        self.LFPM = [None for i in range(self.m)]


    def Add(self, val, t):
        # adds to the shll structure. assuming valid input
        # t: the time of the recieved packet/word val (tk)
        # h: D -> {0,1} ^ 64 bits  (instead of using 32 bits in paper!)
        # x = h(v)
        # i = <x1 x2 ... xb> (first b bits of x), index of 'bucket'/register
        # w = <xb+1 xb+2 ... >
        # <t, p(w)>  (tk, p(wk))
        
        # sha1 takes *bytes* input. converting int,float to str first, so we need to convert str->bytes.
        if isinstance(val,int):
            # convert int to str. can also use an alternate way! (convert into bytes directly)
            # val = val.to_bytes(*size* = 32/64/.. , 'big') 
            val = str(val)
        if isinstance(val,float):
            #converts float to str. gives us up to first 15 decimals after point.
            # self question: can we convert to bytes as an alternate? even if yes, no need because we are hashing the int/float anyway.
            val = str(val)
        if isinstance(val, str):
            val = val.encode('utf-8')


        # hashing the input word - here we used sha1 - can use sha256 too!
        #x = struct.unpack('!Q', sha1(val).digest()[:8])[0] - needs to import struct
        # we use this line now - we can use 'little' instead of big endian. (in practice we take the 64 upper bits of the hash!)
        # note: hex(x) gives us x in hex. hex length = 4 bit so we can make sure that hex(x) is of length 16
        x = int.from_bytes(sha1(val).digest()[:8],'big')
        # obtaining i: by bitwise and - stripping the first b bits only! (m-1 = 2^b-1 = 11..1 for b bits!)
        i = x & (self.m - 1)
        # w - is the rest of the word. shifting the word b bits to the right.
        w = x >> self.b
        # calculating p(w)=Rk, paying attention that p(w) is b bits shorter than 64 because we truncated
        p_w = calculate_p_w(w, 64 - self.b)


        tmp = [] # the updated LFPM list
        
        if self.LFPM[i] is not None:
            for ti, R in self.LFPM[i]:
                # discard old packets
                if (ti < t - self.W):
                    continue
                # discard packets with smaller R
                #note: the below "continue" can be replaced by break. that is because all next R values in
                # the list, are strictly smaller than this R which is <= p_w, therefore we can skip to end and append t,p_w
                if (R <= p_w):
                    continue
                
                tmp.append((ti,R))
        
        tmp.append((t,p_w))
        self.LFPM[i] = list(tmp)
        # note we user list() to create a copy of the list. so won't get deleted after we exit the method

    
    def Merge(self, shll_2):
        # Merges this shll object with another shll_2 object.(Optional)
        # this shll is updated only, while shll_2 is not.
        if type(shll_2) != SlidingHyperLogLog:
            raise TypeError("Cannot merge shll_2 since it is not a Sliding HyperLogLog Object")
        if self.m != shll_2.m:
            raise ValueError("Two Sliding HyperLogLog Objects should have the same number of registers")
            
        for i in range(self.m):
            # nothing to merge (shll 2 is empty)
            if shll_2.LFPM[i] is None:
                continue
            # self is empty, shll 2 is not. assign directly
            if self.LFPM[i] is None:
                self.LFPM[i] = list(shll_2.LFPM[i])
                continue
            # here we have to merge both of them.
            Rmax = None
            tmax = None
            tmp = [] # new list to be built
            merged = list(heapq.merge( *( [self.LFPM[i]] + [shll_2.LFPM[i]] ) ))
            # in order for merge to work, we send *iterables thus the *
            # the reason for coating of each list [lfpm[i]] instead of lfpm[i]:
            # because if we use operator + without outer [], it appends the lists then sends to the function. [tmp1 elements, tmp2 elements]
            # usage of [tmp] + [tmp2] appends in the following way: [ [tmp inner elements] , [tmp2 elements]  ]
            
            for t,R in reversed(merged):
                if tmax is None:
                    tmax = t
                if (t < tmax - self.W):
                    break
                if Rmax is None or R > Rmax:
                    Rmax = R
                    tmp.append((t,R))
                
                # note unlike add/card functions where it was enough to break at first 'success' (suitable R value)
                # here we have to iterate all the list. For example: shll_1 R's: (in order) 6,5,4, shll_2 R's: 5,4,3 
                # after merging 6,5,5,4,4,3 (if we assume suitable overlapping t's) note that since we merge by t's we can have also 5,6,4,3,5,4 etc.
                # in any merge, the inner order of the sublist merged (t's order) is preserved! (e.g. shll_2 cannot be mapped for example 5,..,3,..,4)
                #because we iterate on reversed, we must continue to the end of the list to insert all appropriate R's 
                #because older R's must be STRICTLY increasing (see alternative Add function - similar note)
                    
            tmp.reverse()
            self.LFPM[i] = list(tmp) if tmp else None
        
     
    
    def calculate_cardinality_buckets(self, M):
        # helper function
        # Z_inv = calculate the INVERSE of the indicator Z (Z definition can be shown in paper)
        # Z = 1 / Z_inv, which is what we do when computing E (E = alpha_m * m * m * Z)
        Z_inv = sum(math.pow(2.0, -x) for x in M)
        E = self.alpha_m * float(self.m ** 2) / Z_inv
        
        if (E <= 2.5*self.m):
            # small range correlation
            # count number or registers equal to 0
            V = M.count(0)
            if V > 0:
                E = self.m * math.log(self.m / float(V))
        
        # Did not add large range correction from paper. Redundant because of the use of 64-bit hash functions. 
        return round(E)
  
        
    def EstimateCardinality(self, t, w = 0):
        # t is current timestamp, w is the window (last w units of time)
        # assume t is valid input.
        # if wrong w arg. is sent
        if not 0 < w <= self.W:
            w = self.W
        
        # M is a register-like array of length m
        M = [0 for i in range(self.m)]
        
        # for each lfpm, calculate highest R among the appropriate packets.
        for i in range(self.m):
            if self.LFPM[i] is None:
                continue
            Rmax = 0
            for ti, R in self.LFPM[i]:
                if (ti >= t-w) and R > Rmax:
                    Rmax = R
                    # note: after setting Rmax=R, we can type break into this line, and preserve correctness.
                    # explained: R's in the list are strictly decreasing. therefore, the FIRST R 
                    # that satisfies the time window equation (ti >= t-w), is the representitive R we need
                    # because subsequent R's are strictly smaller!
            M[i] = Rmax
        
        E = self.calculate_cardinality_buckets(M)
        return E


    def EstimateCardinality_list(self, t, w_list):
        # t is current timestamp, w_list contains w's (last w units of time)
        # note: function sorts w_list ascending and returns output accordingly!
        # assume t is valid input.
        if len(w_list) == 0:
            raise ValueError("w_list should not be empty")
        
        # if wrong w arg. is sent
        for i in range(len(w_list)):
            if not 0 < w_list[i] <= self.W:
                w_list[i] = self.W
        #remove duplicates and sort ascending
        w_list = sorted(set(w_list))
        
        n = len(w_list) #for comfort
        
        w_buckets = [None for i in range(n)]
        for i in range(len(w_buckets)):
            # M is a register-like array of length m
            w_buckets[i] = [0 for j in range(self.m)]
 
  
        # for each lfpm, calculate highest R among the appropriate packets. in a manner for all w's sent. explained further in notes
        for k in range(self.m):
            if self.LFPM[k] is None:
                continue
            
            Rmax = 0
            i = 0
            
            for ti, R in reversed(self.LFPM[k]):
                while i <= n-1:
                    if ti >= t-w_list[i]:
                        break # seeking the correct w_i for the w_list that we are in right now. ti (got from the lfpm) is a bad name, don't get confused.                   
                    w_buckets[i][k] = Rmax # filling the array so there won't be middle cells that stay 0 when passing
                    i = i+1
                    # Rmax above can be different tha zero. e.g w_list 1000,100,10, lfpm entries 500,5 so 10-100 has a gap.
                    # w_buckets[i] is a bucket list (M), for time w_i.
                    # we insert in position k because we are now in k-th lfpm.
                
                if i == n:
                    break # out of boundries
                
                if R > Rmax:
                    #side note: if we get here, condition will always be true (since reversed R's are strictly ascending)
                    Rmax = R
                    w_buckets[i][k] = Rmax
                
            while (i <= n-1):
                #refill remaining w's if list is over, but didn't reach all w's
                w_buckets[i][k] = Rmax
                i = i+1
            # the comp. case, finishing all lfpm but entries are very new compared to w's, we have zeros count

        w_results = [0 for i in range(n)]
        for i in range(n):
            w_results[i] = self.calculate_cardinality_buckets(w_buckets[i])
        return w_results


