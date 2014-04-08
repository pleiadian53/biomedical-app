'''
Created on Nov 22, 2013

@author: barnett

This module computes optimal bipartite match in terms 
of either minimum cost or maximum profit. 

Bipartite matching algorithms are given in the module: 
maxWBiMatch
'''

# choose the max-weight bipartite matching algorithm
from maxWBiMatch import maxProfitMatching, \
        minCostMatching, evalMatch
        
# data processor
from preprocess import process_data

### global variables
# [note] 1. If global variables are not preferable, they  
#        can be easily translated into static class 
#        variables by turning affinity module into a class
#        definition. But logically speaking, static class 
#        variables are still "global" within the class definition.  
Vowels = set(['a', 'e', 'i', 'o', 'u'])
DrugSet, ProteinSet = ([], []) # [1] 
         
def countChar(name):
    """
    Given a string (e.g. drug name), count the number of vowels,  
    consonants and other characters (numbers, hyphens, etc). 
    
    Return a summary in 3-tuple: (nVowel, nConsonant, nOther)
    
    [note] 1. not necessary if drug name is guaranteed to have 
              no special chars 
    """
    import re
    # filter non-alphabet characters [1] 
    _name = re.sub('[^a-zA-Z]', '', name)
    nAlphabets, nTotal = (len(_name), len(name))
    
    nVowels = 0
    for ch in _name: 
        if ch.lower() in Vowels: 
            nVowels += 1 
            
    return (nVowels, nAlphabets-nVowels, nTotal-nAlphabets)    
        
     
# compute weight matrix according to the binding rules
def evalWeights(_pdata=True):
    """
    Evaluate weight matrix as an input for a given 
    max weight bipartite matching algorithm such as 
    Hungarian algorithm. 
    
    [note] 1. unnecessary if # of drugs == # of proteins
           2. even length rule: 
                BA = # of vowels * 2
           3. odd length rule: 
                BA = # of consonants * 2.5
           4. increase BA by 25% if any common 
              factors found
           5. or use gcd from fractions module
    """
    global ProteinSet, DrugSet
    
    def gcd(a, b): # [5] 
        # Calculate the Greatest Common Divisor of a and b.
        while b:
            a, b = b, a%b
        return a
    
    # process input data and cache them for later use
    ProteinSet, DrugSet = process_data()
    
    N = max(len(ProteinSet), len(DrugSet))  # [1] 
    w = [[0 for _ in range(N)] for _ in range(N)] 

    for i, protein in enumerate(ProteinSet): 
        for j, drug in enumerate(DrugSet):
            np, nd = (len(protein), len(drug))
            if np % 2 == 0:   # even 
                w[i][j] = countChar(drug)[0] * 2  # [2]
            else:  # odd 
                w[i][j] = countChar(drug)[1] * 2.5  # [3]
            
            if gcd(np, nd) > 1: # [4]
                w[i][j] *= 1.25
    return w

# [test]
def evalRandomAssignment(W=None, _debug=0):
    """
    Evaluate sum of affinity values with 
    random assignment strategy. 
    """
    import random
    
    if W is None: W = evalWeights()
    N = len(W)
    
    random.seed()
    candidates = range(0, N)
    Mu, Mv = ({}, {})
    for i in range(N):
        m = random.choice(candidates)
        Mu[i] = m
        Mv[m] = i
        candidates.remove(m)
    assert len(Mu) == len(Mv), "[Match] Not a perfect match."
    if _debug: 
        for i, j in Mu.items():
            print "protein %d: %s -> drug %d: %s | ba=%f" % \
               (i, ProteinSet[i], j, DrugSet[j], W[i][j])
    
    return (_format(Mu), evalMatch(Mu, W))

def evalOptAssignment(W=None, match_func=maxProfitMatching, _flip=False):
    """
    Evaluate optimal bipartite matching given the weight matrix W. 
    
    *match_func: a bipartite matching function that returns 
                 matching result in 3-tuple: 
                 (Mu, Mv, value)
        where Mu is a dictionary mapping from any two 
                 indep sets of objects e.g. proteins to drugs 
              Mv is a inverse mapping e.g. drugs to proteins 
              value is the sum of all matched weights 
                   e.g. sum of affinity values
    """
    if W is None: W=evalWeights()
    if not hasattr(match_func, '__call__'): 
        raise ValueError, "[evalOptAssignment] Invalid match function: %s" % \
                   str(match_func) 
    Mu, Mv, val = match_func(W, _flip=_flip)
    return (_format(Mu), val)
    
def _format(M):
    """
    Convert matching result from dictionary to a list of tuples.
    """
    assert type(M) == type({})
    return [(i, j) for i, j in M.items()]
        
        
def timeMatching(W=None, match_func=maxProfitMatching, _flip=False):
    """
    Time the given matching algorithm, maxProfitMatching by default. 
    
    *match_func: a function that implements max weight bipartite 
                 matching algorithm. 
                   e.g. maxProfitMatching, minCostMatching, etc. 
                 See module maxWBiMatch
    *_flip: if True, convert the weight matrix to either 
               go from profit to cost matrix or 
               cost to profit matrix. 
    """
    from timer import Timer
    if W is None: W = evalWeights()
    
    Mu, Mv, val = ({}, {}, 0)
    if not hasattr(match_func, '__call__'): 
        raise ValueError, "[timeMatching] Invalid match function: %s" % str(match_func)
    
    try:
        with Timer() as t:
            Mu, Mv, val = match_func(W, _flip=_flip) 
    finally:
        print("> %s took %.03f sec." % (match_func.__name__ + '()', t.interval))
        print "  + assignment: %s" % _format(Mu)
        print "  + sum of affinity: %f" % val
    return (_format(Mu), val) 

def testCountChar(_type='drug'):
    process_data()
    if _type.startswith('d'): dataset = DrugSet
    else: dataset = ProteinSet
    for _str in dataset: 
        print "%s -> %s" % (_str, str(countChar(_str)))
    return

def benchmark():
    import numpy as np
    from datastruct import LinkedLoop
    from timer import Timer
    
    W = evalWeights()  
    
    print "1. Time the computation of matching algorithms ...\n"
    Mu, val = timeMatching(W, maxProfitMatching)

    # wrapper over C++ impl should run faster
    # [note] set _flip to True to convert profit to cost
    Mu2, val2 = timeMatching(np.array(W), minCostMatching, _flip=True)    
    
    # check if there exist different assignments between 
    # two implementations; assignments may not be unique
    print "\n2. Compare the matching results from two different implementations ...\n"
    lloop = LinkedLoop(10)
    for i, pair in enumerate(Mu):
        if Mu2[i][1] != pair[1]: 
            lloop.push("(%d -> %d) vs (%d -> %d)" % \
                          (pair[0], pair[1], Mu2[i][0], Mu2[i][1]))
    if LinkedLoop.total:
        print "  + Examples of different assignment:"
        for e in lloop.content():
            print "     ++ %s" % str(e)    
           
    print "\n3. Now, compare optimum and random assignments ...\n"
    nTrials = 100
    avals = []
    for i in range(nTrials): 
        Mu3, val3 = evalRandomAssignment(W)
        avals.append(val3)
    print "  + random assignment:" 
    print "    ++ pairs: %s" % Mu3
    
    # import math
    # myu = sum(avals)/(len(avals)+0.0)
    # print math.sqrt(sum([pow(e-myu, 2) for e in avals])/(len(avals)+0.0))
    
    print "    ++ max: %f, min: %f, avg: %f, std: %f" % \
         (max(avals), min(avals), sum(avals)/len(avals), np.std(avals))
    
    #print "> history:\n%s\n" % avals
    
    Mu4, val4 = evalOptAssignment(W)
    print "  +  opt assignment:"
    print "     ++ pairs: %s" % Mu4
    print "     ++ value: %f" % val4
     
    

if __name__ == "__main__":
    #testCountChar()
    # test_process_data()
    benchmark()
