import time
import sys
import math, re

def progress_bar(NUM_PROGRESS, NUM_PROGRESS_BIN=0.2, STR_SYS_SYMBOL="=", STR_DES="Progress"):
    NUM_SYM = int(NUM_PROGRESS / NUM_PROGRESS_BIN)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:20s}\>]{1:4.1f}% {2:s}'.format(STR_SYS_SYMBOL*NUM_SYM, NUM_PROGRESS*100, STR_DES))
    sys.stdout.flush()

def loop_progress_cal(ARR_INDEX, ARR_INDEX_MAX, NUM_CUM_MAX=1, NUM_CUM_IND=1, NUM_TOTAL_MAX=1):
    """ Please list from smallest to largest, i.e.: x->y->z """
    if len(ARR_INDEX) == len(ARR_INDEX_MAX):
        for i, i_index in enumerate(ARR_INDEX):
            NUM_IND_PER = (i_index+1)/float(ARR_INDEX_MAX[i])
            NUM_TOTAL_MAX = NUM_TOTAL_MAX * ARR_INDEX_MAX[i] 
            if i >0: NUM_CUM_MAX = NUM_CUM_MAX * ARR_INDEX_MAX[i-1]
            NUM_CUM_IND = NUM_CUM_IND + NUM_CUM_MAX * i_index 
        return NUM_CUM_IND / float(NUM_TOTAL_MAX)
    else:
        print("Wrong dimenstion for in put ARR_INDEX ({0:d}) and ARR_INDEX_MAX ({1:d})".format(len(ARR_INDEX), len(ARR_INDEX_MAX)))

