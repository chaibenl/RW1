#!/usr/bin/python
#used to calculate the read coverage of gene expression
import sys, string

def parseHits(fH):
	line = fH.readline()
	while line:
		try:
			cols = line.strip().split('\t')
			ref = cols[2].strip().split('.')[0]#the contig name
			rstart = int(cols[3])
			rstop = rstart + 34
			read = cols[4]
			occ = int(cols[6]) #number of occurrences the same read is mapped to other place
			try:
				errors = len(cols[7].split(','))#number of mismatches
			except IndexError:
				errors = 0
			if (occ > 0) or (read.count('N') >= 1) or errors > 1: #aligned to multiple genome locations or contain mmore than one 'N's
				low[ref] += 1 #another low quality read
				line = fH.readline()
			else:	
				N[ref] += 1#quality reads mapped to the genome
				q = 1
				feats = ftH[ref]#the genes in this contig
				featList = feats.keys()	
				for ID in featList:
					feat = feats[ID][:2]
					if (rstart >= min(feat) and rstart <= max(feat)):# or (rstop >= min(feat) and rstop <= max(feat)):#within the feature range
						ftH[ref][ID][4] += 1 #add one hit to the feature hash	
						mapped = 1 #this read mapped to at least one feature
						
				line = fH.readline()
		except IndexError:
			line = fH.readline()
	return 1

def readFT(file):#read coordinates of the features
	lines = file.readlines()
	ftHash = {}#feature hash
	for line in lines[1:]:
		cols = line.strip().split('\t')
		ID = cols[0].strip()
		start = int(cols[1])
		end = int(cols[2])
		length = abs(end - start + 1)
		geneName = cols[7]
		contig = cols[-1].split(':')[-1].strip()
		if not contig in ftHash.keys():
			ftHash[contig] = {}
		ftHash[contig][ID] = [start, end, '%s'%length, geneName, 0]

	return ftHash



#Main#
f1 = open(sys.argv[1], 'r')#gene feature table in tab-delimited format
ftH = readFT(f1)
f2 = open(sys.argv[2], 'r')#RNA-seq reads mapped to the genome in Bowtie format
expr = sys.argv[3]+ '_count.txt'#file name to hold gene counts of mapped reads to each feature
GEI = sys.argv[3] + '_summary.txt' #file name to hold gene expression index

N = {}#total number of reads mapped to genomes
GF = {}#reads mapped to gene features
low = {}#reads that did not pass quality filter: unique alignment location and number of "N"s
for contig in ftH.keys():
	N[contig] = 0#total mappable reads
	#GF[contig] = 0 #total reads that are mapped to gene features
	low[contig] = 0#low quality reads
rs = parseHits(f2)#parse the bowtie alignment

out1 = open(expr, 'w')#file to hold the annotated gene expression data
out2 = open(GEI, 'w')

refs = ftH.keys()

out1.write('Contig\tGeneID\tProduct\tReads\n')
out2.write('Contig\tQualityReads2Genome\tQualityReads2Gene\tQualityIntergenic\tLowQuality\n')
for ref in refs:
	featCount = 0 #total reads mapped to features
	fIDs = ftH[ref].keys()
	for fID in fIDs:#iterate feature IDs
		geneLen, name, count = ftH[ref][fID][2:]
		featCount += count
		outString = string.join([ref, fID, name, geneLen, '%s'%count], '\t')
		out1.write(outString + '\n')
	n = N[ref]
	c = featCount
	inter = n - c #intergenic reads: read quality genome mapped minus reads feature mapped	
	l = low[ref]#low quality reads
	outStr2 = string.join([ref, '%s'%n, '%s'%c, '%s'%inter, '%s'%l], '\t')
	out2.write(outStr2 + '\n')

out1.close()
out2.close()
