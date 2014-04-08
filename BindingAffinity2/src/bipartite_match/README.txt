
Introductions to bipartite_match package

@author: Barnett Chiu

[Usage] 
   
   <command> 
       1) Basic usage
       
          python match.py <protein-file> <drug-file>
   
              where <protein-file> holds newline-delimited 
                                   protein names 
                    <drug-file> holds new-delimited drug names 
                    
       2) Set the argument 'versose' to True in main() 
          of match module, 
          i.e. main(verbose=True), to run the benchmark 
          version of the application. This will show
          running times and other details including the
          comparison in different implementations of
          matching algorithms, and between the optimal 
          matching algorithm and random assignments, etc. 
         
   <platform> 
   
      1. This package was tested on the following platform:
      
         a. 2.4 GHz Intel Core i7 (with 8GB memory)
         b. Mac OS X (Version 10.7.5) 
         c. Python 2.7 
            Note that for older versions of Python (prior to 
            Python 2.5), the with statement used in 
            the Timer module may not be supported.
            
      2. external library used: 
         a. numpy
         
         ... and optionally
         b. hungarian 0.2.3 (see module maxWBiMatch or section I)
         
   <optional supportive files> 
      
      1. proteins.txt
      2. drugs.txt 
      

I. Introduction 
   
   The protein-drug matching problem is an instance of maximum-weight 
   bipartite matching, which is known to be solved in polynomial time 
   using Kuhn-Munkres algorithm (aka Hungarian algorithm). One can also 
   formulate it using integer linear programming optimization, or ILP. 
   
   This package attempts to demonstrate and compare two different 
   implementations of Hungarian algorithm with one of them directly 
   utilizing the existing library from pypi.python.org (see 
   the module maxWBiMatch for more details). 
   
   The package consists of the following modules: 
      1) match: the main entry point for this application. 
      2) maxWBiMatch: contains implementations of Kuhn-Munkres algorithm. 
      3) affinity: uses maxWBiMatch to solve max-weight bipartite 
                   matching between proteins and drugs; 
                     
                   A weight matrix W is used as an input to matching 
                   functions. Each weight encodes a binding affinity
                   value, which is determined via the three
                   matching rules that are also defined in this module. 
                   
      4) hungarian.so: a wrapper around a C++ implementation of    
                       Hungarian algorithm compiled from the following
                       package hungarian 0.2.3 from PyPi: 
                       
                       http://pypi.python.org/pypi/hungarian
                       
                       This module is used for benchmarking and 
                       comparison purposes.
                       
      5) preprocess: takes on two input files, one for protein names 
                     and one for drug names, and parses and converts
                     them into a list of protein names and drug names 
                     of string types. 
      
   Other supportive files: 
      6) datastruct.py
      7) timer.py: computes approximate running time for the matching 
                   algorithms. 
                 
   I.1 Weight Matrix:    
   
       There two flavors of weight matrix, W: cost and profit. 
       If the weight represents cost, then we define the optimal 
       assignment in terms of minimum cost. Conversely, if 
       the weight represents profit, then the optimal assignment 
       is defined via maximizing profit. 
       
       The protein-drug matching domain fits with maximum-profit 
       model. However cost matrix and profit matrix are complements 
       to each other and therefore can be easily converted back and
       forth. See flip() function defined in maxWBiMatch for an 
       example. 
       
       hungarian.so uses minimum-cost weight matrix and
       as a result, using it to solve protein-drug assignment requires
       "flipping" the weight (or affinity) matrix (i.e. from profit to 
       cost).  
       
II. ILP
   
    The test cases in maxWBiMatch::demo() are based on the calculation 
    using the procedure given here.   
   
    Max weight bipartite matching can be posed as an instance of ILP: 
    
       minimize      sum<i,j> ( c<i,j> * x<i,j> ) 
       subject to    sum<i> x<i,j> = 1  # one drug to one protein only
                     sum<j> x<i,j> = 1  # one protein to one drug only
                     x<i,j> in {0, 1}
                  
    where c<i,j> represents a cost value
    
    Again, the cost matrix in the above formulation can be 
    converted to profit matrix that we need for making protein-drug 
    assignments, in which case we will be maximizing 
    the objective function instead. For instance, one 
    could perform the following conversion: 
    
       1) find the maximum entry in the matrix 
       2) subtract every cost entry from the maximum
       
    Actually, using W' = -W will suffice since minimum cost is 
    equivalent to maximum profit in terms of the assignment.             
                  
    The test example given in module maxWBiMatch is computed by 
    solving this ILP via the following algorithmic steps: 
    
    e.g. Suppose that we already have weight matrix W available.
    
         W = [[11, 7, 10, 17, 10], 
              [13, 21, 7, 11, 13],
              [13, 13, 15, 13, 14],
              [18, 10, 13, 16, 14], 
              [12, 8, 16, 19, 10]] 
              
         step 1: row and column "reduction" 
                i.e. 
                subtract each row the minimum of each row 
                and subtract each column the minimum of each column 
                to get the following reduced matrix: 
                
                W1 = [[4, 0,  3, 10, 2], 
                      [6, 14, 0, 4,  5], 
                      [0, 0,  2, 0,  0], 
                      [8, 0,  3, 6,  3], 
                      [4, 0,  8, 11, 1]] 
                      
                Test cases: using W1 and W should result in the same 
                            assignment strategy. 
                            (see demo() in maxWBiMatch)
                       
         LOOP through Step 2 ~ Step 4 until full assignments are found 
                      
         step 2: make "attempt" allocations
               
               W1' = [[4, *,  3, 10, 2], 
                      [6, 14, *, 4,  5], 
                      [*, x,  2, x,  x], 
                      [8, x,  3, 6,  3], 
                      [4, x,  8, 11, 1]] 
                     
               from W1, one could make following assignments: 
                  [(0, 1), (1, 2), (2, 0)]
                  
               i.e. if a given row or column has only 1 zero
                    then make the corresponding assignment (denoted
                    in '*'); after the assignments are made, cross out
                    the infeasible entry of 0s (denoted by 'x')
                    
                    in this case, we can only made 3 assignments.
                    therefore, we still need to continue with Step 3 
                    until we find 5/full assignments ...  
                    
          step 3: find minimum unassigned weight: epsilon
          
                3.1) check unassigned rows (in W1, this would be row 
                     3 and 4, the last two rows) 
                  
                loop until no more checking is possible: 
                   3.2) if a checked row has a zero, then check the 
                        corresponding column
                   
                    3.3) if a checked column has a zero, then check the
                        corresponding row 
                        
                3.4) draw "imaginary" lines through unchecked rows and 
                     checked columns
                     
                running 3.1 ~ 3.4 with W1, column 1 (zero-based index) and
                row 1, 2 will have imaginary lines through. All these lines 
                cover all the 0s in W1. We find that the minimum remaining 
                entry, epsilon, is 1
                
          step 4: reduce/transform the weight matrix further with epsilon
                  found in step 3. 
                  
                  foreach entry in W: 
                     subtract epsilon from the entry if No line through
                     add epsilon to the entry if TWO lines through 
                  
                  after the loop above, we get: 
                
                  W2 = [[3, 0,  2,  9,  1], 
                        [6, 15, 0,  4,  5], 
                        [0, 1,  0,  0,  0], 
                        [7, 0,  2,  5,  2], 
                        [3, 0,  7, 10,  0]]
                        
          Notice that W2, when compared to W1, one "new" zero entry 
          is created (not necessarily one extra zero entry though).
                        
          Now we go back to step 2 and make an attempt assignment
          again ... 
                  
                  W2' = [[3, *,  2,  9,  1], 
                         [6, 15, *,  4,  5], 
                         [*, 1,  2,  x,  x], 
                         [7, x,  2,  5,  2], 
                         [3, x,  7, 10,  *]]
                         
          Now, we have 4 assignments: 
                  
              [(0, 1), (1, 2), (2, 0), (4, 4)]
                  
          But still missing 1 more needed assignment
                  
          At this point, W2' by Step 3, has 4 imaginary 
          lines and epsilon = 1
                  
          With epsilon=1, we get another reduced matrix: 
                  
                  W3 = [[2, 0,  1,  8,  0], 
                        [6, 16, 0,  4,  5], 
                        [0, 2,  2,  0,  0], 
                        [6, 0,  1,  4,  1], 
                        [3, 1,  7, 10,  0]]
                        
          W3 again gives us only 4 assignments; thus 
          we need to continue to loop through step 2
          to step 4 until 5 assignments are found
            
          Finally, we will get 
            
           
          W_final = [[0, 0,  0,  6,  0], 
                     [5, 17, 0,  3,  5], 
                     [0, 4,  3,  0,  0], 
                     [4, 0,  0,  2,  1], 
                     [1, 1,  6,  8,  0]]  
          and 
          
          W_final' = [[*, x,  x,  6,  x], 
                      [5, 17, *,  3,  5], 
                      [x, 4,  3,  *,  x], 
                      [4, *,  0,  2,  1], 
                      [1, 1,  6,  8,  *]]
                      
          W_final' gives us full assignments: 
          
          [(0, 0), (1, 2), (2, 3), (3, 1), (4, 4)], 
          
          with which we then obtain the optimal cost
          by referring to initial W: 
          
          11 + 7 + 13 + 10 + 10 = 51
          
          This is answer is justified by 
          demo() in maxWBiMatch. 
          
          
          
          
       
             
                  
                  
                  
                  
                
                
                
                     
                  
                     
                    
                    
                    
               
                       
               
    
    
    
    
                            

   

         

