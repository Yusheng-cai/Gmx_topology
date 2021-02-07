from topology.sep_top import *
import sys

if __name__ == "__main__":
    fname = sys.argv[1] 
    t = topology(fname)

    t.write_ff("test.itp")

    
