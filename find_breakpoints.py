#!/usr/bin/env python

# requires 2 arguments:
#   string path to sam file
#   name of output bed file you want to make (will append '_hits.bed' and '_rejects.bed' to the file name)
#
# example: python find_soft_clipped_ends.py 'my_alignments.sam' 'my_softclips'



import re
import operator
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Find the soft-clipped breakpoints')
    
    parser.add_argument('input', type=str, help='SAM format alignment file')
    parser.add_argument('output', type=str, help='Name for output file(s), _hits.bed will be appended')
    parser.add_argument('-t','--threshold', type=int, required=False, default=3, help='Threshold for including hits in the hits file. Default = 3.')
    parser.add_argument('-w','--window', type=int, required=False, default=3, help='Sliding window size for averaging number of hits at a position. Default = 5.')
    
    return parser.parse_args()

def find_clip_locations():
    args = parse_args()
    if not args.window % 2 == 1:
        print 'Window size must be an odd number.'
        return
    
    # Create an empty breakpoint dictionary
    # This will store all soft-clipped breakpoint hits
    # It will be a dictionary of dictionaries - one for each chromosome:
    ## {chr1 : {position1 : #hits, position2 : #hits} , chr2 : {postion1 : #hits, position2 : #hits}}
    breakpoint_counts = {}

    with open(args.input, 'r') as in_file:
        for line in in_file:
            # split fields
            entries = line.split('\t')
            # check if it's a header line, and if so, skip it
            if re.search('^@[A-Z][A-Z]$', entries[0]):
                pass
            else:
                #grab the chromosome name, position, and the cigar string
                chr_name, read_start, cigar = entries[2], int(entries[3]), entries[5]
                #parse the cigar info, exclude hard-clipped regions (these can only be on the very edges, outside of soft-clips if soft-clipping is present)
                map_regions = re.findall('[0-9]+[MIDNSP=X]', cigar)

                # Get the first and last items from the cigar string and see if it's soft-clipped (last letter = S). Soft clips will only ever be at the ends (inside would be I/D/N)
                # If so, find the breakpoint and add it to the dictionary
                if map_regions[0][-1] == 'S':
                    breakpoint = read_start
                    try:
                        breakpoint_counts[chr_name][breakpoint] = breakpoint_counts[chr_name].get(breakpoint, 0) + 1
                    ## if the dictionary doesn't exist yet for this chromosome, make it, then add the breakpoint
                    except KeyError:
                        breakpoint_counts[chr_name] = {}
                        breakpoint_counts[chr_name][breakpoint] = breakpoint_counts[chr_name].get(breakpoint, 0) + 1
                if map_regions[-1][-1] == 'S':
                    breakpoint = read_start + len(entries[9]) - int(map_regions[-1][:-1])
                    try:
                        breakpoint_counts[chr_name][breakpoint] = breakpoint_counts[chr_name].get(breakpoint, 0) + 1
                    ## if the dictionary doesn't exist yet for this chromosome, make it, then add the breakpoint
                    except KeyError:
                        breakpoint_counts[chr_name] = {}
                        breakpoint_counts[chr_name][breakpoint] = breakpoint_counts[chr_name].get(breakpoint, 0) + 1                
                            

    # Iterate over the keys and find the sliding-window average for each position.  Keep the breakpoint if it's >= threshold.
    # First grab the list of chromosome dictionaries, then for each chromosome, iterate over all stored positions in that chromosome
    with open(args.output + '_hits.bed', 'w') as out_hits, open(args.output + '_rejects.bed', 'w') as out_rejects:
        for chr_name in breakpoint_counts.keys():
            hits_list = []
            rejects_list = []
            chr_dict = breakpoint_counts[chr_name]
            for position in chr_dict.keys():
                # Add up all the positions +/- half of the sliding window width, then average them, to determine whether or not to keep the hit
                average_hits = sum([chr_dict.get((position + i - args.window/2), 0) for i in range(args.window)]) / float(args.window)
                # Add the position and the number of hits to the appropriate list
                if average_hits >= args.threshold:
                    hits_list.append( (position, chr_dict[position]))
                else:
                    rejects_list.append( (position, chr_dict[position]) )

            # sort by position of breakpoint and write to files
            hits_list.sort(key=operator.itemgetter(0))
            rejects_list.sort(key=operator.itemgetter(0))
    
            for i in hits_list:
                out_hits.write(chr_name + '\t' + str(i[0]) + '\t' + str(i[0] + 1) + '\t' + 'softclip\t' + str(i[1]) + '\n')
            for i in rejects_list:
                out_rejects.write(chr_name + '\t' + str(i[0]) + '\t' + str(i[0] + 1) + '\t' + 'softclip\t' + str(i[1]) + '\n')

            
if __name__ == '__main__':
    find_clip_locations()
