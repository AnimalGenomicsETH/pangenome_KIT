#!/usr/bin/env python3

import bisect
import sys

def translate_bed_to_graph(gfa,bed):
    
    boundaries = [int(line.rstrip().split()[5][5:]) for line in open(gfa,'r') if line[0] == 'S' and 'SR:i:0' in line]

    with open(bed,'r') as fin:
        for line in fin:
            s_coord, e_coord = (int(i) for i in line.rstrip().split()[1:3])
            gene = line.rstrip().split()[3]

            s_node = bisect.bisect_left(boundaries,s_coord)
            s_offset = s_coord - boundaries[s_node-1]
            
            e_node = bisect.bisect_left(boundaries,e_coord)
            e_offset = e_coord - boundaries[e_node-1]
            
            if s_node == e_node:
                print(f's{s_node}\t{s_offset}\t{e_offset}\t{gene}')
            else:
                for i in range(s_node,e_node+1):
                    if i == s_node:
                        print(f's{s_node}\t{s_offset}\t{boundaries[s_node]-boundaries[s_node-1]}\t{gene}')
                    elif i != e_node:
                        print(f's{i}\t{0}\t{boundaries[i]-boundaries[i-1]}\t{gene}')
                    else:
                        print(f's{e_node}\t{0}\t{e_offset}\t{gene}')

if len(sys.argv) == 3:
    translate_bed_to_graph(*sys.argv[1:3])
else:
    print('python translate_bed_to_graph.py <gfa> <bed>\n !! Incorrect arguments !!\nPrints output to stdout')
