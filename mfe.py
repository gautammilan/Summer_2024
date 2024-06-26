from collections import defaultdict
import math
from linearparition import  VALID, get_value, beamprune, seq

R= 8.3144
T= 1

# print(len(seq))
def mfe(seq,b):

    n= len(seq)
    Q= [defaultdict(lambda: float("inf")) for _ in range(n+1)]
    back= defaultdict(list)
    for j in range(1, n+1): Q[j-1][j]= 0
    for j in range(1, n+1):
        for i in Q[j-1]:
            if Q[j][i]> Q[j-1][i]:
                Q[j][i]= Q[j-1][i]
                back[i,j]= [(i,j-1)]

            pair, value= get_value(i-1,j, seq)
            if pair: 
                for k in Q[i-2]:
                    temp= Q[i-2][k]+Q[j-1][i]+value
                    if Q[j][k]> temp:
                        Q[j][k]=temp
                        back[k,j]= [(k,i-2), (i,j-1)]
        beamprune(Q,j, b)
    
    def solution(h_edges):
        if len(h_edges)==1:return solution(back[h_edges[0]])+"."
        elif len(h_edges)==2: return solution(back[h_edges[0]])+"("+solution(back[h_edges[1]])+ ")"
        else: return ""
        
    return Q[n][1], solution(back[1,n])
# seq= "GAUGCCGUGUAGUCCAAAGACUUC"
print(len(seq))
print("GEH ",mfe(seq, 100))