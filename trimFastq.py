#!/usr/bin/python
#This script trims off 3' end of Illuminia reads by Read Segment Quality Control indicator on standard Phred+33 scale: ASCII 33-71 (in which the Read Segment Quality Control indicator is designated as "#", while it is "B" on Phred+66 scale)
import sys, string
if len(sys.argv) != 3:
	print 'Usage: trimFastq.py seq.fastq minLength[INT]'
	sys.exit()

f = open(sys.argv[1], 'r').read()#the fastq file
min = int(sys.argv[2])#the minimum length required
hist = {} #count of each length
recs = f.split('@HWI')[1:]
root = sys.argv[1].split('.')[0]
output1 = open(root + '_trimmed.fastq', 'w')
for rec in recs:
	lines = rec.split('\n')
	ID = '@HWI' + lines[0].strip()#sequence ID
	seq = lines[1].strip()#sequence string
	sep = lines[2].strip()#separator
	qual = lines[3].strip()#quality characters
	
	length = len(seq)
	poss = range(length)
	poss.reverse()
	cut = length
	for pos in poss:
		Q = qual[pos]
		if Q == '#':
			cut = pos	
		else:
			break
	if not hist.has_key(cut):
		hist[cut] = 1
	else:
		hist[cut] += 1
	if cut  < min:
		continue

	out = [ID, seq[:cut], sep, qual[:cut]]
	output1.write(string.join(out, '\n') + '\n')

output1.close()
keys = hist.keys()
keys.sort()
output2 = open(root + '_trimmed.sum', 'w')
output2.write('Length\tCount\n')
for length in keys:
	out = ['%s'%length, '%s'%hist[length]]
	output2.write(string.join(out, '\t') + '\n')
output2.close()

