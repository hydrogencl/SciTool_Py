import time
import sys
import math, re

def progress_bar(NUM_PROGRESS, NUM_PROGRESS_BIN=0.2, STR_SYS_SYMBOL="=", STR_DES="Progress"):
    NUM_SYM = int(NUM_PROGRESS / NUM_PROGRESS_BIN)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:20s}\>]{1:4.1f}% {2:s}'.format(STR_SYS_SYMBOL*NUM_SYM, NUM_PROGRESS*100, STR_DES))
    sys.stdout.flush()
    
