import math
from collections import defaultdict
import random
from linearparition import inside, VALID, get_value, R, T


T= 0.1
def prod(lst):
    val= 1
    for v in lst:
        val= val*v
    return val

def binary_search_first_greater(lst, target):
    if not lst:
        return None

    low = 0
    high = len(lst) - 1
    result = None

    while low <= high:
        mid = (low + high) // 2
        if lst[mid][1] > target:
            result = lst[mid]  
            high = mid - 1 
            if high >= 0 and lst[high][1] < target:
                return result[0]  
        else:
            low = mid + 1 
    return result[0]  


class LinearSampling:
    def __init__(self, seq,k):
        self.seq= seq
        self.k= k
        self.S= defaultdict(list)
    def get_weight(self,subs):  
        if len(subs)==1:
            return 0
        else:
            _, s2= subs
            return get_value(s2[0]-1,s2[1]+1,seq)[1]

    def in_edges(self,span):
        i,j= span
        edges= []
        if j== i-1: return [((i,j-1), [])]
        if i in self.Q[j-1]: edges.append(((i,j),[(i,j-1)]))
        for k in self.Q[j-1]:
            if k-1>= i and i in self.Q[k-2]:
                if self.seq[k-1-1]+self.seq[j-1] in VALID:
                     edges.append(((i,j), [(i,k-2), (k,j-1)]))
        return edges
    
    def sum_edges(self,v, saving):
        s = 0
        for e in self.in_edges(v):
            # print(e)
            _,subs= e
        
            if subs:
                energy= self.get_weight(subs)
                s += math.exp((energy/(R*T) * prod([self.Q[j][i] for i,j in subs])))
                if saving: self.S[v].append((e,s))
        return s

    def sample_lazy(self,v):
        i,j= v
        if v not in self.visited:
            self.sum_edges(v, True)
            self.visited.add(v)
       
        h_edges= self.S[v]
        if h_edges:
            print(h_edges, self.Q[j][i])
            t= random.uniform(0, self.Q[j][i])
            print(h_edges, self.Q[j][i],t)
            _,subs= binary_search_first_greater(h_edges,t)
            # print("the selected subs", subs)
            if len(subs)==1:return self.sample_lazy((i,j-1))+"."
            elif len(subs)==2: return self.sample_lazy(subs[0])+"("+self.sample_lazy(subs[1])+ ")"
        else: return ""
    

    def main_lazy(self):
        samples= []
        self.Q= inside(self.seq, self.k)
        self.visited= set()
        n = len(self.seq)

        for _ in range(self.k):
            samples.append(self.sample_lazy((1, n)))
        return samples

seq= "GCACG"
LS= LinearSampling(seq, 6)
lr= LS.main_lazy()
print(lr)