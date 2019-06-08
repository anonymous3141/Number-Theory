#euclid, extend euclid, and chinese remainder
import random,math
array=[]
e=2.7182818284590452353602874713
pi=3.141592653589793238462643383
I=1j #constant definitions
def gcd(a,b,arrayname=False): #a>=b for all a,b and input, boo arrayname=array to dump q values
    if b==0:
        return(a)
    else:
        r=a%b
        if arrayname:
            exec(arrayname+".append(int((a-a%b)/b))")
        return(gcd(b,r,arrayname))
    
def extend_euclid(a,b): #calculates inverse of a in a=1 mod b with b>=a, gives bugged output if there is no inverse e.g. 10 mod 12 has no inverse
    global q
    q=[]    
    p=[0,1]
    gcd(b,a,"q")
    #print(q)
    for c in range(2,len(q)+1):
        p.append((p[c-2]-p[c-1]*q[c-2])%b) #the algorithm
    #print(p)
    if (p[-1]*a)%b==1: #checked if output correct i.e if there exists inverse
        return(p[-1])
    return(False) #no inverse

def chinese_remainder(a,m):    #solves set of modular equations with all modulos coprime
    #a is the set of remainders and m is the set of modulos: i.e all equations is: x=a mod m
    #returns an list of two elements, solution x and modulo M in that order
    M=1
    y=[] #inverses
    for c in m:
        M*=c
    for k in range(len(a)):
        y.append(extend_euclid((M/m[k])%m[k],m[k])) #find all y_k where y_k is the inverse of M/m_k mod m_k
    x=0 #x is the solution (modulo M is computed at end)
    for c in range(len(a)):
        x+=int(a[c]*M/m[c]*y[c]) #compute sigma solution
    return([x%M,M]) #solution is always modulo M if one exists

class polynomial:
    def __init__(self, coeffs):
        self.coeffs=coeffs
        self.degree=len(coeffs)-1
        self.dydx=[]
    def p(self,value):
        out=0
        for c in range(self.degree+1):
            out+=value**(self.degree-c)*self.coeffs[c]
        return(out)
    def ddx(self):
        arr=[]
        for c in range(self.degree):
            arr.append(self.coeffs[c]*(self.degree-c))
        self.dydx=arr[:]
    def pprime(self,value): #deriative of p(x) at a point
        out=0
        for c in range(self.degree):
            out+=value**(self.degree-c-1)*self.dydx[c]
        return(out)
    def findroot(self,x=0):
        while True:
            x=0
            while abs(self.p(x))>10**(-8):
                if self.pprime(x)==0:
                    x+=1
                    break
                x=x-self.p(x)/self.pprime(x)
            return(x)
        
    def findvalue(self,goal):
        a0=None
        b0=None
        for c in range(-10000,10000):
            if self.p(c)>=goal:
                a0=c
            else:
                b0=c
        return(self.bisection(a0,b0,goal))
                
            
    def bisection(self,a,b,goal): #function f
        #print(a,b)
        c=(a+b)/2
        if abs(self.p(c)-goal)<10**-8:
            return(c)
        else:
            #print(abs(self.p(c)-goal),c,self.p(c))
            if self.p(c)>goal: #given a>=b
                return(self.bisection(c,b,goal))
            else:
                return(self.bisection(a,c,goal))
        
def integral(low,up,f):
    area=0
    boo=-1
    if None not in [low,up]:        
        for c in range(500*low,500*up+1):
            area+=f(c/500)*0.002

    else:
        if low != None:
            boo=1
            current=low
        else:
            current=up
        while True:
            a=0
            for c in range(current*500,current*500+500*boo,boo):
                a+=f(c/500)*0.002
            area+=a
            current+=1
            if abs(a)<=10**(-8):
                print(a,"derp")
                break
    return(area)

def g(x):
    return(e**(x))
def modexp(a, b, p): #b is power of 2
    if b==1: return(a)
    else:
        return(modexp(a,b//2,p)**2%p)
def primitive(p, k): #calculate a primitive 2^kth root of unity mod p
    #precondition: 2^k divides p-1
    for c in range(50):
        a=random.randint(1,p-1)
        if modexp(a,(p-1)//2,p) != 1:
            #print(a)
            return(modexp(a,(p-1)//(1<<k),p))
def fft(a,b=True): #if b=False then IFFT is computed instead, without dividing by N
    #recursive, length of a is pow of 2
    if len(a)==1:
        return(a)
    else:
        a0=[a[i] for i in range(0,len(a),2)]
        a1=[a[i] for i in range(1,len(a)+1,2)]
        y0=fft(a0,b)
        y1=fft(a1,b)
        Y=[0]*len(a)
        W=1
        Wn=math.cos(2*pi/len(a))+I*math.sin(2*pi/len(a))
        if not b:
            Wn=1/Wn
        for k in range(len(y0)):
            t=W*y1[k] #butterfly
            Y[k]=y0[k]+t
            Y[k+len(a)//2]=y0[k]-t
            W=W*Wn
        return(Y)

def ifft(y):
    a=fft(y,False)
    for k in range(len(a)): #divide by n
        a[k]=a[k]/len(a)
    return(a)
    
def multiply(P,Q):
    #P,Q are polynomial coefficient vectors of length N=2^n, same length
    #not important as just pad with zeros if needed
    P2=[0]*len(P)
    Q2=[0]*len(Q)
    P2.extend(P)
    Q2.extend(Q) #pad with zeros so that we have 2N values for our fft
    p=fft(P2)
    q=fft(Q2) #point representation
    r=[p[i]*q[i] for i in range(0,2*len(P))] #point interpolate
    #poly[i] is the ith power term
    return(ifft(r))
    
def inverse(a,b,p): #calculate x with ax equiv b mod prime p
    inva=1; e=1
    while e<p-1:
        if e & (p-2):
            inva=(inva*modexp(a,e,p))%p
        e*=2
    return((inva*b)%p)
                
        
            
"""print(gcd(5**63+1,4**63))
print(extend_euclid(1008,2017))
print(chinese_remainder([1,4,6],[3,5,7])) 
print(chinese_remainder([2,3],[3,8]))
#print(bin("herp zerp gerp")>>10)
p=polynomial([3,0,-1,1])
print(p.p(10))
print(p.ddx())
print(p.dydx,p.degree)
print(p.pprime(1))
print(p.findroot())
print(p.findvalue(50))
print(integral(None,0,g))
print(ifft(fft([420,69,666,-69,-666,-420,888,8888])))
print(multiply([1,3,3,1],[1,3,3,1]))
print(inverse(420,69,29))"""
