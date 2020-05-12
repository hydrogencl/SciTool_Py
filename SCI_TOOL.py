import time
import sys
import math, re

def progress_bar(NUM_PROGRESS, NUM_PROGRESS_BIN=0.2, STR_SYS_SYMBOL="=", STR_DES="Progress"):
    NUM_SYM = int(NUM_PROGRESS / NUM_PROGRESS_BIN)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:20s}\>]{1:4.1f}% {2:s}'.format(STR_SYS_SYMBOL*NUM_SYM, NUM_PROGRESS*100, STR_DES))
    sys.stdout.flush()
    
def loop_progress_cal(ARR_INDEX, ARR_INDEX_MAX, NUM_OUT=1, NUM_ERROR=-9999.999):
    if len(ARR_INDEX) == len(ARR_INDEX_MAX):
        for i, i_index in enumerate(ARR_INDEX):
            NUM_OUT = NUM_OUT * (i_index/float(ARR_INDEX_MAX[i])) 
        return NUM_OUT
    else:
        print("Wrong dimenstion for in put ARR_INDEX ({0:d}) and ARR_INDEX_MAX ({1:d})".format(len(ARR_INDEX), len(ARR_INDEX_MAX)))
        return NUM_ERROR
