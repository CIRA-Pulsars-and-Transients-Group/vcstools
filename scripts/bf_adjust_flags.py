#!/usr/bin/env python

msg='''
This is a little script that will adjust the flags.txt as output by get_delays.
 It takes the flagged_tiles.txt (the name is not important as long as there is one tile per line
 tile numbers are assumed to be zero based.) as created by the one who ran the RTS to find a 
 calibration solution. Based on this file it will adjust the flags in flags.txt (as output by get_delays,
 one entry per line, two entries per tile) as needed. 
'''
def modify_flags(rts_flag_file, flags_file):
    try:
        with open(rts_flag_file,'r') as rts_flags_file:
            flagged_tiles = [int(flag.strip()) for flag in rts_flags_file.readlines()]
    except:
        print "Couldn't read flagged_tiles."
        return 1
    try:
        with open(flags_file,'r') as flags_f:
            flags = flags_f.readlines()
    except:
        print "Couldn't read flags.txt"
        return 1
    try:
        # tile numbering is zero based (hence can be used as list index) and
        # there is an entry for each tile and polarisation in flags.txt from get_delays.
        for tile in flagged_tiles:
            flags[tile*2] = '0.0\n'
            flags[tile*2 + 1] = '0.0\n'
        with open(flags_file,'w') as flags_f:
            for flag in flags:
                flags_f.write(flag)
        return 0
    except:
        print "Couldn't modify flags.txt"
        return 1


if __name__ == '__main__':
    from sys import argv,exit
    from os import path
    if len(argv) < 2 or argv[1] == '-h':
        print "Usage: {0} /path/to/flagged_tiles.txt /path/to/pointingDir/chan/flags.txt".format(argv[0])
        print msg        
        exit(1)
    rts_flag_file = argv[1]
    flags_file = argv[2]
    if not path.exists(rts_flag_file):
        print "Cannot open {0}. Will leave flags untouched. Aborting.".format(rts_flag_file)
        exit(1)
    if not path.exists(flags_file):
        print "Cannot open {0}. Will leave flags untouched. Aborting.".format(flags_file)
        exit(1)
    exit(modify_flags(rts_flag_file, flags_file))
