from collections import defaultdict
import math
import random

R= 8.3144
T= 0.01
VALID= {"AU", "UA", "CG", "GC", "GU", "UG"}

def get_value(i,j, seq):
    i,j= i-1,j-1

    if seq[i]+seq[j] in ["AU", "UA"]: return True, -1
    elif seq[i]+seq[j] in ["CG", "GC"]: return True, -1
    elif seq[i]+seq[j] in ["GU", "UG"]: return True, -1
    else: return False, None


NEG_INF = float('-inf')

def logsum(x, y):
    if x < y:
        x, y = y, x  # Swapping x and y
    if y > NEG_INF / 2 and x - y < 11.8624794162:
        x = math.log1p(math.exp(x-y)) + y
    return x


def selecttopb(candidates,b):

    def qselect(array,k):
        try:
            if len(array) <= 1:
                return array[0]
        except:
            print(array,k)
            exit()
        pivot_index = random.randint(0, len(array) - 1) 
        pivot = array[pivot_index]

        left = [x for i, x in enumerate(array) if x < pivot and i != pivot_index]
        right = [x for i, x in enumerate(array) if x >= pivot and i != pivot_index]

        k_adj = k - 1

        if k_adj < len(left):
            return qselect(left,k)
        elif k_adj == len(left):
            return pivot
        else:
            return qselect(right,k_adj - len(left))

    value= qselect(candidates.values(), b)
    new_candidates=  {key:val for key, val in candidates.items() if val<value}
    
    if len(new_candidates)<b:
        for key, val in candidates.items():
            if val==value:
                new_candidates[key]= value
            if len(new_candidates) == b:
                break
    return new_candidates

def beamprune(Q, j, b):
    candidates = {}

    # First gather candidates with safe access to Q[i-1][1] if i-1 exists in Q and 1 exists in Q[i-1]

    for i in Q[j]:
        if 1 in Q[i-1]:
            candidates[i] = Q[i-1][1] * Q[j][i]
    # print("inside beamprune")
    if b < len(candidates):
        # print(len(candidates))
        candidates = selecttopb(candidates, b)  
    keys_to_delete = [i for i in Q[j] if i not in candidates]

    for i in keys_to_delete:
        if j==len(seq) and i==0 and i not in keys_to_delete:
            print("DEleting I which is wrong")
        del Q[j][i]

def inside(seq,b):
    n= len(seq)
    Q= [defaultdict(float) for _ in range(n+1)]
   
    for j in range(1, n+1): Q[j-1][j]= 1 
    for j in range(1, n+1):
        for i in Q[j-1]:
            Q[j][i]+=  Q[j-1][i]
            pair, value= get_value(i-1,j, seq)
            if i-2>-1 and pair: 
                for k in Q[i-2]:
                    Q[j][k]+= Q[i-2][k]*Q[j-1][i]*math.exp(-value/(R*T))
        if j!=len(seq):
            beamprune(Q,j, b)
    return Q

def outside(seq, Q):
    n= len(seq)
    print("the value of partition function is", Q[n][1])

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
                    p[i-1,j]+= (Q_out[k,j]*Q[i-2][k]*math.exp(-energy/(R*T))*Q[j-1][i])/Q[n][1]           
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
                    p[i-1,j]= logsum(p[i-1,j], (Q_out[k,j]+Q[i-2][k]+(-energy/(R*T))+Q[j-1][i]-Q[n][1]))
    return p, Q_out

# seq= "GAUGCCGUGUAGUCCAAAGACUUC"
# Q= inside(seq, len(seq))
# p,Q_out= outside(seq, Q)
# print(Q[len(seq)][1])



# Q_= log_inside(seq, 100)
# p_,Q_out_= log_outside(seq, Q_)

# print(Q_[len(seq)][1])
# print(math.exp(Q_[len(seq)][1]))

# print(p)
