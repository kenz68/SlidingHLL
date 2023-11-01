
import math
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




class HyperLogLog:
    """ Implementation of HyperLogLog class.
    Properties of the class:
    m = 2^b - the number of registers
    b -  (log 2 of m above)
    alpha_m - the const used for corrction of hash bias (read article)
    M - an array of m registers, used as in the article
    """
    
    def __init__(self, param):
        """
        The constructor of the class.
        We enable the end-user to use the constructor in one of two ways:
        1. Sending the log of the number of required registers directly: 
        param needs to be of type int, b = log_2(m) with b in range [4..16]
        2. Sending the allowed Std. error as a fraction (0 < r < 1)
        param then is of type float!
        Std. error is defined as in the article: (E' - E) / E where E is actual cardinality, E' is our estimate
        """
        
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

        #common code for both cases: set b, m, alpha_m, init. m registers
        self.b = b
        # m = 2 ** b (2^b)
        self.m = 1 << b 
        self.alpha_m = calculate_alpha_m(b)
        # M(1)... M(m) = 0 // m registers initialized            
        self.M = [ 0 for i in range(self.m) ]



    def Add(self, val):
        # adds to the hll structure. assuming valid input
        # h: D -> {0,1} ^ 64 bits  (instead of using 32 bits in paper!)
        # x = h(v)
        # i = <x1 x2 ... xb> (first b bits of x), index of 'bucket'/register
        # w = <xb+1 xb+2 ... >
        # M[i] = max(M[i], p(w))
        
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
        # updating the corresponding register, paying attention that p(w) is b bits shorter.
        self.M[i] = max(self.M[i], calculate_p_w(w, 64 - self.b))
        
    
    def Merge(self, hll_2):
        # Merges this hll object with another hll_2 object.(Optional)
        # this hll is updated only, while hll_2 is not.

	# in practice we might need to merge n counters. to write another funciton is redundant.
	# maxing (1,2,...,n) reg vals is of same time complexity as maxing (1,2) into (1) then (1,3) , ... , (1,n)
	# because in each case we have n comparisons


        if type(hll_2) != HyperLogLog:
            raise TypeError("Cannot merge hll_2 since it is not a HyperLogLog Object")
        if self.m != hll_2.m:
            raise ValueError("Two HyperLogLog Objects should have the same number of registers")
        self.M = [max(self.M[i],hll_2.M[i]) for i in range(self.m) ]

    
    def EstimateCardinality(self):
        # Z_inv = calculate the INVERSE of the indicator Z (Z definition can be shown in paper)
        # Z = 1 / Z_inv, which is what we do when computing E (E = alpha_m * m * m * Z)
        Z_inv = sum(math.pow(2.0, -x) for x in self.M)
        E = self.alpha_m * float(self.m ** 2) / Z_inv
        
        if (E <= 2.5*self.m):
            # small range correlation
            # count number or registers equal to 0
            V = self.M.count(0)
            if V > 0:
                E = self.m * math.log(self.m / float(V))
        
        # Did not add large range correction from paper. Redundant because of the use of 64-bit hash functions. 
        return round(E)
        
        