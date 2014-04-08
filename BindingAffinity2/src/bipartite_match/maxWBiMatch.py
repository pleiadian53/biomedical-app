#!/usr/bin/python

"""
Kuhn-Munkres, aka the Hungarian algorithm. Complexity O(n^3)
Computes a max weight perfect matching in a bipartite graph
for min weight matching, simply negate the weights.

@author: barnett 
@reference: 1. ecole polytechnique - c.durr - 2009
            2. https://pypi.python.org/pypi/hungarian
            3. https://pypi.python.org/pypi/munkres


Global variables:
       n = number of vertices on each side
       U,V vertex sets
       lu,lv are the labels of U and V resp.
       the matching is encoded as 
       - a mapping Mu from U to V, 
       - and Mv from V to U.
    
    The algorithm repeatedly builds an alternating tree, rooted in a
    free vertex u0. S is the set of vertices in U covered by the tree.
    For every vertex v, T[v] is the parent in the tree and Mv[v] the
    child.

    The algorithm maintains minSlack, s.t. for every vertex v not in
    T, minSlack[v]=(val,u1), where val is the minimum slack
    lu[u] + lv[v] - w[u][v] over u in S, and u1 is the vertex that
    realizes this minimum.

    Complexity is O(n^3), because there are n iterations in
    maxWeightMatching, and each call to augment costs O(n^2). This is
    because augment() makes at most n iterations itself, and each
    updating of minSlack costs O(n).
"""

import numpy as np
import sys

def improveLabels(val):
    """ 
    Change the labels, and maintain minSlack. 
    """
    for u in S:
        lu[u] -= val
    for v in V:
        if v in T:
            lv[v] += val
        else:
            minSlack[v][0] -= val

def improveMatching(v):
    """ 
    Apply the alternating path from v to the root in the tree. 
    """
    u = T[v]
    if u in Mu:
        improveMatching(Mu[u])
    Mu[u] = v
    Mv[v] = u

def slack(u,v): return lu[u]+lv[v]-w[u][v]

def augment():
    """ 
    Augment the matching, possibly improving the labels on the way.
    """
    while True:
        # select edge (u,v) with u in S, v not in T and min slack
        ((val, u), v) = min([(minSlack[v], v) for v in V if v not in T])
        assert u in S
        if val>0:        
            improveLabels(val)
        # now we are sure that (u,v) is saturated
        assert slack(u,v)==0
        T[v] = u                            # add (u,v) to the tree
        if v in Mv:
            u1 = Mv[v]                      # matched edge, 
            assert not u1 in S
            S[u1] = True                    # ... add endpoint to tree 
            for v in V:                     # maintain minSlack
                if not v in T and minSlack[v][0] > slack(u1,v):
                    minSlack[v] = [slack(u1,v), u1]
        else:
            improveMatching(v)              # v is a free vertex
            return

def maxProfitMatching(weights, _flip=False):  # minimum cost
    """ 
    Compute best assignment of maximum profit; i.e. each weight 
    represents profile (rather than cost).  
    
    Given w, the weight matrix of a complete bipartite graph,
    returns the mappings Mu : U->V ,Mv : V->U encoding the matching
    as well as the value of it.
    
    *_flip: if True, convert input weight matrix to cost matrix 
    """
    global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
    
    w  = np.array(weights)
    if _flip: 
        w = -w   # more profit is equivalent to less cost
        
    n  = len(w)
    U  = V = range(n)
    lu = [ max([w[u][v] for v in V]) for u in U]  # start with trivial labels
    lv = [ 0                         for v in V]
    Mu = {}                                       # start with empty matching
    Mv = {}
    while len(Mu)<n:
        free = [u for u in V if u not in Mu]      # choose free vertex u0
        u0 = free[0]
        S = {u0: True}                            # grow tree from u0 on
        T = {}
        minSlack = [[slack(u0,v), u0] for v in V]
        augment()
    
    # val. of matching is total edge weight 
    val = sum(lu)+sum(lv)
    
    return (Mu, Mv, 0.0-val if _flip else val)

def minCostMatching(weights, _flip=False):
    """
    Compute minimum-cost assignment (by default) where 
    weights represents cost by default.

    *_flip: if True, convert cost matrix *weights 
               to an equivalent profit matrix
    """
    _library= False
    try: 
        import hungarian
        _library = True
    except: 
        import os
        msg = "[minCostMatching] hungarian.so not found in %s?" % \
               os.getcwd()
        print msg
    
    _weights = np.array(weights)
    if not _library: 
        if not _flip: # W represents cost by default 
            _weights = flip(weights) # convert to profit
        else: # W would represent profit had we wanted to convert it
            pass
        return maxProfitMatching(_weights)
    
    # library hungarian is available, use it
    if _flip:  # convert to max profit problem
        _weights = flip(weights)
    
    match1, match2 = hungarian.lap(_weights)
    
    Mu, Mv = ({}, {})
    for i, j in enumerate(match1):
        Mu[i] = j
    for i, j in enumerate(match2):
        Mv[j] = i
        
    # evaluate total cost (or profit) using the original weight matrix
    # assert evalMatch(Mu, weights) == evalMatch(Mv, weights), "Inconsistent total weight"
    return (Mu, Mv, evalMatch(Mu, weights))
    
def evalMatch(M, W, _T=False):
    _M = dict(M)
    _W = np.array(W)
    #print "dim: %s" % str(_W.shape)
    if _T: _W = _W.transpose()
    totalWeight = 0
    for i, j in _M.items():
        totalWeight += _W[i][j]    
    return totalWeight

def flip(W):
    """
    Convert a cost matrix an equivalent profit matrix and 
    vice versa; i.e. if W is a cost matrix then output is 
    profit matrix. 
    """
    weights = np.array(W)
    upperbound = np.max(weights)
    _weights = []
    for row in weights:
        cost_row = []
        for col in row:
            #cost_row += [sys.maxsize - col]
            cost_row += [upperbound - col]
        _weights += [cost_row]
    return _weights
  
def demo():
    #import numpy as np
 
    #W = [[1,2,3,4],[2,4,6,8],[3,6,9,12],[4,8,12,16]] # profits, max profit = 30
    #W =  [[5,9,3,6],[8,7,8,2],[6,10,12,7],[3,10,8,6]]  # costs, min-cost = 18
    W = [[11, 7, 10, 17, 10], 
         [13, 21, 7, 11, 13],
         [13, 13, 15, 13, 14],
         [18, 10, 13, 16, 14], 
         [12, 8, 16, 19, 10]]     # costs, min-cost = 51
    
    _flip = False
    feedback = raw_input('> W represents cost? (y/n)  \n')
    if feedback.lower().startswith('y'): 
        _flip = True
    
    W = np.array(W)
    Mu, Mv, val = maxProfitMatching(W, _flip=_flip)  # convert to min. cost
    print "> Using maximum profit match ...\n"
    print "(%s, %s, %f)" % (Mu, Mv, val)
    
    # adding or subtracting a common weight from 
    # rows or columns should not change the assignment
    print "\n> Subtracting a common weight from rows and columns ...\n"
    N = len(W)
    _W = np.array(W, copy=True)
    for i in range(N):
        mval = min(_W[i,:])
        for j in range(N): 
            _W[i,j] -= mval   # row "reduction"
    Mu2, Mv2, _ =  maxProfitMatching(_W, _flip=_flip)
    
    print "(%s, %s, %f)" % (Mu2, Mv2, evalMatch(Mv2, W, _T=True))
    
    Wt = _W.transpose()
    for i in range(N):
        mval = min(Wt[i,:])
        for j in range(N): 
            Wt[i,j] -= mval   # column "reduction"
    _W = Wt.transpose()
    Mu3, Mv3, _ = maxProfitMatching(_W, _flip=_flip)
    print "(%s, %s, %f)" % (Mu3, Mv3, evalMatch(Mv3, W, _T=True))
    
    for match in [Mu2, Mu3, ]:
        for i, j in Mu.items():
            assert match[i] == j, "Inconsistent matching after reduction!"

    print "\n> Compare the result with the impl from PyPi (which uses min-cost match) ...\n"
    
    Mu4, Mv4, val = minCostMatching(W, _flip=(not _flip))
    print "(%s, %s, %f)" % (Mu4, Mv4, val)
    
    return 

  
if __name__ == "__main__":
    demo()