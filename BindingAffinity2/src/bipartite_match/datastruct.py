'''
Created on Nov 22, 2013

This module implements LinkedLoop representing 
a circular queue. 

@author: barnett
'''

class LinkedLoop(object):
    head = 0
    tail = 0
    total = 0
    def __init__(self, capacity=10):
        self.size = capacity
        self._list = [None] * capacity
        
    def push(self, item):
        self._list[LinkedLoop.tail] = item
           
        LinkedLoop.tail = (LinkedLoop.tail+1) % self.size 
        if self.isFull(): 
            LinkedLoop.head = LinkedLoop.tail
          
        if LinkedLoop.total < self.size: 
            LinkedLoop.total += 1 
        return                         
    def enqueue(self, item):
        return self.append(item)
    
    def isEmpty(self):
        return LinkedLoop.total == 0
    def isFull(self):
        return LinkedLoop.total == self.size
    
    def pop(self):  # FIFO
        """
        Pop the element in FIFO manner
        """
        if LinkedLoop.total == 0: 
            return None
            # raise RuntimeError, "Empty!"
             
        item = self._list[LinkedLoop.head]
        self._list[LinkedLoop.head] = None
        LinkedLoop.head += 1
        if LinkedLoop.head == self.size: 
            LinkedLoop.head = 0
            
        LinkedLoop.total = LinkedLoop.total-1
        return item
    def dequeue(self):
        return self.pop()
    
    def content(self):
        return self._list
    def __str__(self):
        return str(self._list) 
       
if __name__ == "__main__":
    N = 5 
    lloop = LinkedLoop(N)
    
    M= 27
    print ">>>> cycle 1 <<<<"
    for i in range(M):
        lloop.push(i)
        print lloop, lloop.total
    for j in range(N): 
        a = lloop.pop()
        print "popped %s | content: %s | total: %d" % (a, lloop, lloop.total)
        
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    print "> empty? %s" % bool(lloop.isEmpty())
        
    print ">>>> cycle 2 <<<<"
    for i in range(4):
        lloop.push(20+i); 
        print lloop
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    
    for j in range(4):
        a = lloop.pop()
        print "popped %s | %s" % (a, lloop)
    
    print "> empty? %s" % bool(lloop.isEmpty())    
    
    print ">>>> cycle 3 <<<<"
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    for i in range(M):
        lloop.push(i+30)
        print lloop, lloop.total
    for j in range(N): 
        a = lloop.pop()
        print "popped %s | content: %s | total: %d" % (a, lloop, lloop.total) 
         
         
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    print ">>>> cycle 4 <<<<"
    for i in range(4):
        lloop.push(20+i); 
        print lloop
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    
    for j in range(4):
        a = lloop.pop()
        print "popped %s | %s" % (a, lloop)
    print "> empty? %s" % bool(lloop.isEmpty())  
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)
    
    print "> try popping empty list ..."
    for i in range(3):
        a = lloop.pop()
        print "popped %s | content: %s | total: %d" % (a, lloop, lloop.total)   
    print "> head, tail pointers should not change"
    print "(head: %s, tail: %s)" % (lloop.head, lloop.tail)    
    
            
            
            
        
    

