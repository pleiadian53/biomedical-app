'''
Created on Nov 22, 2013

@author: barnett

Process the input data; i.e. 
   protein file and drug file 

'''
import sys, os
 
if len(sys.argv) != 3:
    msg = "Usage : python %s protein_file drug_file\n" % os.path.basename(sys.argv[0])
    msg += "        where protein_file holds newline-separated protein names\n"
    msg += "              drug_file holds newline-separated drug names\n"
    sys.stderr.write(msg)
    raise SystemExit(1)

### system/global variables 

PROTEIN_FILE = sys.argv[1]
DRUG_FILE = sys.argv[2]
#PROTEIN_FILE = 'proteins.txt'
#DRUG_FILE = 'drugs.txt'
CURDIR = os.getcwd()

    
def _check_files():
    msg = ''
    st = 0
    if not os.path.exists(PROTEIN_FILE): 
        msg += "[Input] Could not find %s in %s" % (PROTEIN_FILE, CURDIR)
        st += 1
    if not os.path.exists(DRUG_FILE):
        msg += "[Input] Could not find %s in %s" % (DRUG_FILE, CURDIR)
        st += 1
    if st: raise RuntimeError, msg
    return

# parse input file

def process_data():
    """
    Read data from files and store them in lists.
    """
    _check_files()
    return ([line.strip() for line in open(PROTEIN_FILE) if line.strip()], 
            [line.strip() for line in open(DRUG_FILE) if line.strip()]) 
   
def test_process_data(_debug=1):
    proteinSet, drugSet = process_data()
    msg = ''
    if _debug: 
        msg += "> # of proteins:\n%d\n" % len(proteinSet) 
        msg += "> proteins:\n%s\n" % proteinSet 
        msg += "> # of drugs:\n%d\n" % len(drugSet) 
        msg += "> drugs:\n%s\n" % drugSet 
        print msg
    
    # assuming that # of drugs and # of proteins are the same
    assert len(drugSet) == len(proteinSet), \
        "[Input] number of proteins %d, number of drugs %s" % (len(drugSet), len(proteinSet))
    return

if __name__ == "__main__":
    test_process_data()