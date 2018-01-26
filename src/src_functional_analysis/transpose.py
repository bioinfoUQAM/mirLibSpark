#!/usr/bin/python
programme = "transpose.py"
author = "Chao-Jung (Julie) Wu"
version = "1.00.00"
date = "2017-01"
'''
This program ROW-TO-COLUMN transposes a tab delimited text file


in.txt
1   2   3   4   5
6   7   8   9   10
11  12  13  14  15


out.txt
1	6	11	
2	7	12	
3	8	13	
4	9	14	
5	10	15	

'''

def transpose_txt(infile, outfile):
    with open(infile, 'r') as f:
        lis = [x.rstrip('\n').split('\t') for x in f]
        #lis = [x.split() for x in f]

    fho = open (outfile, 'w')
    for x in zip(*lis):
        for y in x:
            print(y+'\t', end='', file=fho, flush=True)
        print('', file=fho)
		
