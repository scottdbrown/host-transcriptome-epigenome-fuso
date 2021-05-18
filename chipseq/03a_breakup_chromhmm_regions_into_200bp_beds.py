TITLE = "Breakup ChromHMM regions"
DESC = "Given ChromHMM regions across genome, break up into the 200bp regions and create bed file."
'''
Used as input to GREAT to weight each region by size.

Date: March 17, 2021
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time


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
    parser.add_argument("chromhmm_file", help = "ChromHMM regions file", type = str)
    parser.add_argument("output_file_base", help = "Output ChromHMM regions file (200bp regions) base name", type = str)
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



    ## want matrix of what state each region is:
    #chrom_set_1 = ["chr{}".format(x) for x in range(1,9)]

    to_write = [[]]
    #to_write.append("\t".join(states))
    write_i = 0
    cur_num_reg = 1
    max_num_reg = 400000

    for line in open(args.chromhmm_file, "r"):
        line = line.rstrip().split("\t")
        chrom = line[0]
        start_pos = int(line[1])
        region_length = int(line[2]) - int(line[1])
        num_regions = int(region_length / 200)
        state = line[3]
        named_state = "{}-{}_{}".format(state, chrom, start_pos)
        
        if cur_num_reg >= max_num_reg:
            ## increment write_i
            write_i += 1
            to_write.append([]) ## add empty list
            cur_num_reg = 1

        for i in range(num_regions):
            to_write[write_i].append("{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t65,65,65".format(chrom,
                    start_pos + 200*i,
                    start_pos + 200*(i+1),
                    "{}-{}".format(named_state, i),
                    start_pos + 200*i,
                    start_pos + 200*(i+1)))
            cur_num_reg += 1
        
    
    for i in range(len(to_write)):
        statprint("Writing to {}".format("{}_{}.bed".format(args.output_file_base, i+1)))
        out = open("{}_{}.bed".format(args.output_file_base, i+1), "w")
        out.write("\n".join(to_write[i]))
        out.write("\n")
        out.close()

    statprint("Done.")