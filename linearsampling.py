import math
from collections import defaultdict
import random
from linearparition import inside,outside, VALID, get_value, R, T
import argparse


def prod(lst):
    val= 1
    for v in lst:
        val= val*v
    return val

def binary_search_first_greater(lst, target):
    if not lst:
        return None
    # print(lst,target)
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
    def __init__(self, seq,k, b):
        self.seq= seq
        self.k= k
        self.b= b
        self.S= defaultdict(list)
    def get_weight(self,subs):  
        if len(subs)==1:
            return 0
        else:
            _, s2= subs
            return get_value(s2[0]-1,s2[1]+1,self.seq)[1]

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
                s += math.exp(-(energy/(R*T))) * prod([self.Q[j][i] for i,j in subs])
                if saving: self.S[v].append((e,s))
        return s

    def sample_lazy(self,v):
        i,j= v
        if v not in self.visited:
            self.sum_edges(v, True)
            self.visited.add(v)
       
        h_edges= self.S[v]
        # print(v, h_edges)
        
        if h_edges:
            # print(h_edges, self.Q[j][i])
            t= random.uniform(0, self.Q[j][i])
            # print(h_edges, self.Q[j][i],t)
            _,subs= binary_search_first_greater(h_edges,t)
            # print("the selected subs", subs)
            if len(subs)==1:return self.sample_lazy((i,j-1))+"."
            elif len(subs)==2: return self.sample_lazy(subs[0])+"("+self.sample_lazy(subs[1])+ ")"
        else: return ""
    

    def main_lazy(self, Q):
        samples= []
        self.Q= Q
        self.visited= set()
        n = len(self.seq)

        for _ in range(self.k):
            samples.append(self.sample_lazy((1, n)))
        return samples

def calculate_all_base_pair_probabilities(structures):
    n = len(structures[0])
    total = len(structures)
    
    # Initialize a matrix to store base pair probabilities
    base_pair_probabilities = [[0 for _ in range(n)] for _ in range(n)]
    
    for structure in structures:
        for i in range(n):
            for j in range(i+1, n):
                if structure[i] == '(' and structure[j] == ')':
                    base_pair_probabilities[i][j] += 1
    
    # Convert counts to probabilities
    for i in range(n):
        for j in range(i+1, n):
            base_pair_probabilities[i][j] /= total
    
    return base_pair_probabilities

def parse_structure(structure):
    stack = []
    pairs = {}
    for index, char in enumerate(structure):
        if char == '(':
            stack.append(index)
        elif char == ')':
            if stack:
                opening_index = stack.pop()
                pairs[opening_index] = index
    return pairs


def aggregate_pairs(structures):
    pair_counts = {}
    for structure in structures:
        pairs = parse_structure(structure)
        for i, j in pairs.items():
            if (i, j) not in pair_counts:
                pair_counts[(i, j)] = 0
            pair_counts[(i, j)] += 1
    return pair_counts

def calculate_probabilities(pair_counts, total_structures):
    probabilities = {}
    for pair, count in pair_counts.items():
        probabilities[pair[0]+1, pair[1]+1] = float(count) / total_structures
    return probabilities



def main():
   
    parser = argparse.ArgumentParser(description="Calculate RNA folding probabilities.")
    parser.add_argument("--sequence", type=str, help="The RNA sequence.")
    parser.add_argument("--beam_size", type=int, help="Beam size for the algorithm.")
    parser.add_argument("--sample_size", type=int, help="Number of sample")

    args = parser.parse_args()


    seq= args.sequence
    n= len(seq)

    Q = inside(seq, args.beam_size)
    p, Q_out = outside(seq, Q)

    ls= LinearSampling(seq, args.sample_size, args.beam_size)
    samples= ls.main_lazy(Q)
    pair_count= aggregate_pairs(samples)
    sample_p= calculate_probabilities(pair_count, len(samples))
    print("\n")
    print("VERIFICATION")
    print("\n")
    print("The probability from the parititon function:")
    print("\n")
    print(dict(p))

    print("\n")
    print("The probability from sampling:")
    print("\n")
    print(sample_p)

if __name__ == "__main__":
    main()