import random
import math

#To generate random prime less than N
def randPrime(N):
    primes = []
    for q in range(2,N+1):
        if(isPrime(q)):
            primes.append(q)
    return primes[random.randint(0,len(primes)-1)]

# To check if a number is prime
def isPrime(q):
    if(q > 1):
        for i in range(2, int(math.sqrt(q)) + 1):
            if (q % i == 0):
                return False
        return True
    else:
        return False

#pattern matching
def randPatternMatch(eps,p,x):
    N = findN(eps,len(p))
    q = randPrime(N)
    return modPatternMatch(q,p,x)

#pattern matching with wildcard
def randPatternMatchWildcard(eps,p,x):
    N = findN(eps,len(p))
    q = randPrime(N)
    return modPatternMatchWildcard(q,p,x)

# return appropriate N that satisfies the error bounds
#Proof that the given value of N satisfies error bounds is in the following drive link: https://drive.google.com/file/d/1M1_svIgBTf9QQMp55uXVxdRH1ShMZa1p/view?usp=drivesdk
def findN(eps,m):
    N=100*((m/eps)**2)
    return int(N//1)

# Return sorted list of starting indices where p matches x
#Using horner's method to evaluate f(p) in O(m). Returns f(p) mod q to satisfy logq space constraints.
def horner(st,x,q):
    m=len(st)
    result = (ord(st[0])-ord("A")) 
    for i in range(1, m):
        result = result*x + (ord(st[i])-ord("A"))
    return result%q
#Storing c, which is power of 26 corresponding to first position of pattern or substring, separately(mod q) to avoid computing it n times. fp is f(p), and fj stores the hash function for document substring,
#starting from index 0, all three computations take logq bits of
#space and happen in O(m) time. Then we loop over all possible index values in the document and evaluate fj for all of them, and try to match fj=fp.
#To avoid computing fj separately for all indices, we use f(starting from index j+1)=26*(f(j)-contribution of index getting deleted)+contri from index getting added.
#For example, if document string is "ABCDE" and we need to match substrings of length 4, we'll first compute f("ABCD") and then use f("BCDE")=26(f("ABCD")-26**(3)*0)+4)
#to compute f("BCDE"). This is true because of the polynomial nature of the hash function and ensures that f(j+1) can be computed using f(j)
#in O(log q) time (since all quantities have been taken mod q).
#This means that the for loop executes in O(nlogq) time and the initial computations of fp and fj took O(mlogq) time and O(logq) space.
#This ensures O((m+n)logq) time complexity for the complete function and O(k+logn+logq) space complexity.
def modPatternMatch(q,p,x):
    L=[]
    c=(26**(len(p)-1))%q
    fp=horner(p,26,q)
    fj=horner(x[0:len(p)],26,q)
    for j in range(1,len(x)-len(p)+1):
        if fj==fp:
            L.append(j-1)
        fj=(26*((fj)-(c*(ord(x[j-1])-ord('A'))))+ord(x[j+len(p)-1])-ord('A'))%q
    if fj==fp:
        L.append(len(x)-len(p))
    return L
# Return sorted list of starting indices where p matches x
#Checks if a character is the wildcard character.
def wildcard(ch):
    if ord(ch)<ord("A") or ord(ch)>ord("Z"):
        return True
    return False
#A variation of horner's method which deletes the contribution from the wildcard character by making its contribution 0. For example, "CD2B" will have
#the hash function (26^3)*2+(26^2)*3+(26^1)*0+(26^0)*1. Returns the position of the wildcard character too.
def new_horner(st,x,q):
    m=len(st)
    if wildcard(st[0]):
        po=0
        result=0
    else:
        result = (ord(st[0])-ord("A")) 
    for i in range(1, m):
        if not wildcard(st[i]):
            result = result*x + (ord(st[i])-ord("A"))
        else:
            po=i
            result=result*x
    return (result%q , po)
#Another variation of horner's method which does the same thing as new_horner but takes the position of wildcard as an argument. This is being
#applied for initial computation of the hash function for document string.
def pos_horner(st,x,q,pos):
    m=len(st)
    if pos!=0:
        result = (ord(st[0])-ord("A"))
    else:
        result=0
    for i in range(1, m):
        if pos!=i:
            result = result*x + (ord(st[i])-ord("A"))
        else:
            result=result*x
    return result%q
#A modification of modPatternMatch. Every computation of fj will ignore the contribution from that index of the document string which matches with the index
#at which wildcard is stored. k represents the power of 26 corresponding to wildcard position. The meaning
#of c is the same as modPatternMatch. k and c have been stores separately mod q to avoid repeated computations.
#At each step, first we add contribution of character at wildcard position. Then, we change fj just like in modPatternMatch to represent the next substring, and then
#we delete contribution of character at wildcard position for this new substring.
#For example, if document string is "XYTAG" and we are comparing subtrings of length 4, and wildcard is at index 2 in pattern,
#first we take f("XYAA") as fj (deleting contribution from T. Note that no character is actualy changed in the string). Then we change "XYAA" to "XYTA",
#then, similar to modPatternMatch, change it to "YTAG", and finally delete contribution of character at wildcard index to get f("YTAG").
#This ensures that when we match fp mod q= fj mod q, the wildcard index is ignored
#for both fp and fj and all other indices are matched.
#Note that complexities will remain the same here as modPatternMatch, and since q is O(log(m/eps)), we have time complexity O((m+n)log(m/eps))
#and space complexity O(k+logn+log(m/eps)).
def modPatternMatchWildcard(q,p,x):
    L=[]
    c=(26**(len(p)-1))%q
    y=new_horner(p,26,q)
    fp=y[0]
    pos=y[1]
    k=(26**(len(p)-pos-1))%q
    fj=pos_horner(x[0:len(p)],26,q,pos)
    for j in range(1,len(x)-len(p)+1):
        if fj==fp:
            L.append(j-1)
        fj=(fj+k*(ord(x[pos+j-1])-ord("A")))%q   #Adding contribution of old term at wildcard index.
        fj=(26*((fj)-(c*(ord(x[j-1])-ord('A'))))+ord(x[j+len(p)-1])-ord('A')-(k*(ord(x[pos+j])-ord("A"))))%q #Note the extra subtraction compared to modPatternMatch. That is deleting contribution of wildcard term.
    if fj==fp:
        L.append(len(x)-len(p))
    return L
