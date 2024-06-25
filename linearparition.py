from collections import defaultdict
import math
import random

R= 8.3144
T= float('inf')
VALID= {"AU", "UA", "CG", "GC", "GU", "UG"}

def get_value(i,j, seq):
    i,j= i-1,j-1

    if seq[i]+seq[j] in ["AU", "UA"]: return True, -2
    elif seq[i]+seq[j] in ["CG", "GC"]: return True, -3
    elif seq[i]+seq[j] in ["GU", "UG"]: return True, -1
    else: return False, None

def logsum(a,b):
    try:
        return math.log(math.exp(a)+math.exp(b))
    except:
        return float("-inf")

def selecttopb(candidates,b):

    def quickselect(arr, k):
        if len(arr) == 1:
            return arr[0]
        
        pivot = random.choice(arr)
        
        lows = [el for el in arr if el < pivot]
        highs = [el for el in arr if el > pivot]
        pivots = [el for el in arr if el == pivot]
        
        if k < len(lows):
            return quickselect(lows, k)
        elif k < len(lows) + len(pivots):
            return pivots[0]
        else:
            return quickselect(highs, k - len(lows) - len(pivots))
   
    value= quickselect(candidates.values(), b)
    new_candidates=  {key:val for key, val in candidates.items() if val>value}
    # print(new_candidates)
    while len(new_candidates)<b:
        for key, val in candidates.items():
            if val==value:
                new_candidates[key]= value
    return new_candidates

def beamprune(Q,j,b):
    candidates= dict()
    for i in Q[j]:
        candidates[i]= Q[i-1][1]*Q[j][i]
    # print(j, candidates)
    if b<len(candidates):
        candidates= selecttopb(candidates, b)   
    poss= Q[j].copy()
    for i in poss:
        if i not in candidates:
            del Q[j][i]




def inside(seq,b):
    n= len(seq)
    Q= [defaultdict(float) for _ in range(n+1)]
   
    for j in range(1, n+1): Q[j-1][j]= 1 
    # print("inside ", Q)
    for j in range(1, n+1):
        for i in Q[j-1]:
            Q[j][i]+=  Q[j-1][i]
            pair, value= get_value(i-1,j, seq)
            if i-2>-1 and pair: 
                for k in Q[i-2]:
                    Q[j][k]+= Q[i-2][k]*Q[j-1][i]*math.exp(-value/(R*T))
    beamprune(Q,j, b)
    return Q

def outside(seq, Q):
    n= len(seq)
    Q_out= defaultdict(float)
    p= defaultdict(float)
    Q_out[1,n]= 1
    for j in range(n, 0, -1):
        for i in Q[j-1]:
            Q_out[i, j-1]+= Q_out[i,j]
            pair, energy= get_value(i-1,j,seq) #converting it to zero based index
            if pair and i-2>-1:
                for k in Q[i-2]: 
                    Q_out[k,i-2]+= Q_out[k,j]*Q[j-1][i]*math.exp(-energy/(R*T))
                    Q_out[i,j-1]+= Q_out[k,j]*Q[i-2][k]*math.exp(-energy/(R*T))
                    p[i-1,j]+= (Q_out[k,j]*Q[i-2][k]*math.exp(-energy/(R*T)))/Q[n][1]
    return p, Q_out

def log_inside(seq,b):
    n= len(seq)
    Q= [defaultdict(lambda: float('-inf')) for _ in range(n+1)]
   
    for j in range(1, n+1): Q[j-1][j]= 0 
    for j in range(1, n+1):
        for i in Q[j-1]:
            Q[j][i]= logsum(Q[j][i],Q[j-1][i])
            pair, value= get_value(i-1,j,seq)
            if i-2>-1 and pair: 
                for k in Q[i-2]:
                    Q[j][k]= logsum(Q[j][k], Q[i-2][k]+Q[j-1][i]+(-value/(R*T)))
    beamprune(Q,j, b)
    return Q

def log_outside(seq, Q):
    n= len(seq)
    Q_out= defaultdict(lambda:float('-inf'))
    p= defaultdict(lambda:float('-inf'))
    Q_out[1,n]= 0
    for j in range(n, 0, -1):
        for i in Q[j-1]:
            Q_out[i, j-1]= logsum(Q_out[i, j-1], Q_out[i,j])
            pair, energy= get_value(i-1,j,seq) #converting it to zero based index
            if pair and i-2>-1:
                for k in Q[i-2]: 
                    Q_out[k,i-2]= logsum(Q_out[k,i-2], Q_out[k,j]+Q[j-1][i]+(-energy/(R*T)))
                    Q_out[i,j-1]= logsum(Q_out[i,j-1], Q_out[k,j]+Q[i-2][k]+(-energy/(R*T)))
                    p[i-1,j]= logsum(p[i-1,j], (Q_out[k,j]+Q[i-2][k]+(-energy/(R*T))-Q[n][1]))
    return p, Q_out
seq= "AGGCAUCAAACCCUGCAUGGGAGCG"
Q= inside(seq, len("AGGCAUCAAACCCUGCAUGGGAGCG"))
Q_out= outside(seq, Q)
