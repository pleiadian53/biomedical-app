'''
Created on Nov 22, 2013

@author: barnett 

Main entry for finding optimal assignments 
between proteins and drugs in terms of their 
maximum binding affinity. 

Usage : python match.py file-1 file-2
        where file-1 holds newline separated protein names
              file-2 holds newline separated drug names 
'''

#from affinity import evalOptAssignment, evalRandomAssignment, \
#       timeMatching, evalWeights, benchmark
       
# from maxWBiMatch import maxProfitMatching, minCostMatching
import affinity

def main(verbose=False):
    if not verbose: 
        assignments, value = affinity.evalOptAssignment()
        msg = "> Assignment:\n%s\n" % assignments 
        msg += "> BA value:  \n%f\n" % value
        print msg
    else: 
        affinity.benchmark()
    return

if __name__ == "__main__":
    main()

