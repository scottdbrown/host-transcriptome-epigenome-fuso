TITLE = "Find Consistent Chromatin States"
DESC = "Given ChromHMM states for replicates, identify regions of the genome which have consistent calls"
'''
Date: February 24, 2021
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time
from collections import Counter


## Declare global variables

DEBUG = False
VERB = False

## Classes and functions

class bcolors:
    CYAN = '\033[1;36;40m'
    BLUE = '\033[1;34;40m'
    GREEN = '\033[1;32;40m'
    YELLOW = '\033[1;33;40m'
    RED = '\033[1;31;40m'
    BOLDWHITE = '\033[1;37;40m'
    DARK = '\033[1;30;40m'
    PURPLE = '\033[1;35;40m'
    ENDC = '\033[0m'

def statprint(msg, msg_type = "STATUS"):
    typeColour = ""
    if msg_type == "ERROR":
        typeColour = bcolors.RED
    elif msg_type == "WARNING":
        typeColour = bcolors.YELLOW
    elif msg_type == "DEBUG":
        typeColour = bcolors.GREEN
    elif msg_type == "SUBPROCESS":
        typeColour = bcolors.GREEN
        msg_type = "     " + msg_type
    else:
        typeColour = bcolors.BOLDWHITE

    print("{message_color}{message_type}{end_color} {time_color}[{datetime}]{end_color}: {message}".format(message_color = typeColour, 
             message_type = msg_type,
             end_color = bcolors.ENDC, 
             time_color = bcolors.BLUE, 
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg))


## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--rep1", help = "Replicate 1", type = str)
    parser.add_argument("--rep2", help = "Replicate 2", type = str)
    parser.add_argument("--rep3", help = "Replicate 3", type = str)
    parser.add_argument("--majority_outfile", help = "Bed file with majority state calls", type = str)
    parser.add_argument("--concordant_outfile", help = "Bed file with concordant state calls", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")


    ## Script content here

    ## At a high level, I need to step through the three replicates, line by line, but not necessarily together.
    ## As these are .bed files, each line is a region, which will be sized differently.
    ## So, I need to have 3 open files, and increment the line in each specific file when needed
    ## (when the window is past the region)

    window_size = 200
    window_start = 0
    window_end = window_start + window_size
    window_chr = ""
    this_chr_max = 99999999

    majority_windows = []
    concordant_windows = []

    ## open files
    rep1 = open(args.rep1, "r")
    rep2 = open(args.rep2, "r")
    rep3 = open(args.rep3, "r")

    # move past header
    rep1_line = rep1.readline()
    rep2_line = rep2.readline()
    rep3_line = rep3.readline()

    header_line = rep1_line

    # get first lines
    rep1_line = rep1.readline()
    rep2_line = rep2.readline()
    rep3_line = rep3.readline()
    
    # get first chromosome name from first file.
    window_chr = rep1_line.split("\t")[0]
    statprint("Working on {}".format(window_chr))

    ## while there is data for all 3 samples (because once one has hit the end, we can't do our test)
    while(all([rep1_line, rep2_line, rep3_line])):
        
        ctr = Counter([rep1_line.split("\t")[3], rep2_line.split("\t")[3], rep3_line.split("\t")[3]]).most_common(1)
        ## ctr holds [("state", count)]

        # if all are equal (have concordance)
        if ctr[0][1] == 3:
            concordant_windows.append((window_chr, window_start, window_end, "{}".format(ctr[0][0])))  ## add the tuple
            majority_windows.append((window_chr, window_start, window_end, "{}".format(ctr[0][0])))
        # if majority (>=2) but not concordant (==2)
        elif ctr[0][1] == 2:
            concordant_windows.append((window_chr, window_start, window_end, "{}".format("X_NoConsistentState")))  ## add the tuple
            majority_windows.append((window_chr, window_start, window_end, "{}".format(ctr[0][0])))  ## add the tuple
        else:
            concordant_windows.append((window_chr, window_start, window_end, "{}".format("X_NoConsistentState")))  ## add the tuple
            majority_windows.append((window_chr, window_start, window_end, "{}".format("X_NoConsistentState")))  ## add the tuple
        
        # Increment window
        window_start += window_size
        window_end += window_size

        # check that all regions cover window and are on same chromosome
        # if window is past end of region or on different chromosome
        if window_end > int(rep1_line.split("\t")[2]) or window_chr != rep1_line.split("\t")[0]:
            ## get next region/line
            rep1_line = rep1.readline()
            if rep1_line.split("\t")[0] != window_chr:
                ## the next region caused a shift to the next chromosome
                ## can move the window
                window_chr = rep1_line.split("\t")[0]
                statprint("Working on {}".format(window_chr))
                window_start = 0
                window_end = window_start + window_size
        if window_end > int(rep2_line.split("\t")[2]) or window_chr != rep2_line.split("\t")[0]:
            ## get next region/line
            rep2_line = rep2.readline()
            if rep2_line.split("\t")[0] != window_chr:
                ## the next region caused a shift to the next chromosome
                ## can move the window
                window_chr = rep2_line.split("\t")[0]
                statprint("Working on {}".format(window_chr))
                window_start = 0
                window_end = window_start + window_size
        if window_end > int(rep3_line.split("\t")[2]) or window_chr != rep3_line.split("\t")[0]:
            ## get next region/line
            rep3_line = rep3.readline()
            if rep3_line.split("\t")[0] != window_chr:
                ## the next region caused a shift to the next chromosome
                ## can move the window
                window_chr = rep3_line.split("\t")[0]
                statprint("Working on {}".format(window_chr))
                window_start = 0
                window_end = window_start + window_size
    
    ## reached the end in at least 1 file. Finish up.

    

    ## Parse through stored windows and merge if they are adjacent.
    ## Need to check chromosome, and if current start matches previous end.
    statprint("Starting to write output.")

    ## concordant
    out = open(args.concordant_outfile, "w")
    out.write(header_line)

    prev_chrom = None
    prev_start = None
    prev_end = None
    prev_state = None
    
    for window in concordant_windows:
        chrom, start, end, state = window
        if not prev_chrom:
            ## initialize
            prev_chrom = chrom
            prev_start = start
            prev_end = end
            prev_state = state
        elif chrom == prev_chrom and start == prev_end and state == prev_state:
            ## part of the same region, extend
            prev_end = end
        else:
            ## not part of the same region, write the previous and update.
            out.write("{c}\t{s}\t{e}\t{st}\t0\t.\t{s}\t{e}\t65,65,65\n".format(c=prev_chrom,
            s=prev_start, e=prev_end, st=prev_state))
            
            prev_chrom = chrom
            prev_start = start
            prev_end = end
            prev_state = state
    ## write last region
    out.write("{c}\t{s}\t{e}\t{st}\t0\t.\t{s}\t{e}\t65,65,65\n".format(c=prev_chrom,
    s=prev_start, e=prev_end, st=prev_state))
    
    out.close()

    ## majority
    out = open(args.majority_outfile, "w")
    out.write(header_line)

    prev_chrom = None
    prev_start = None
    prev_end = None
    prev_state = None
    
    for window in majority_windows:
        chrom, start, end, state = window
        if not prev_chrom:
            ## initialize
            prev_chrom = chrom
            prev_start = start
            prev_end = end
            prev_state = state
        elif chrom == prev_chrom and start == prev_end and state == prev_state:
            ## part of the same region, extend
            prev_end = end
        else:
            ## not part of the same region, write the previous and update.
            out.write("{c}\t{s}\t{e}\t{st}\t0\t.\t{s}\t{e}\t65,65,65\n".format(c=prev_chrom,
            s=prev_start, e=prev_end, st=prev_state))
            
            prev_chrom = chrom
            prev_start = start
            prev_end = end
            prev_state = state
    ## write last region
    out.write("{c}\t{s}\t{e}\t{st}\t0\t.\t{s}\t{e}\t65,65,65\n".format(c=prev_chrom,
    s=prev_start, e=prev_end, st=prev_state))
    
    out.close()

    statprint("Done.")