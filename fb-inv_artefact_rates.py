#!/usr/bin/env python
#coding: utf-8


## DEPENDENCIES ##
from cigar import Cigar
import pysam
import argparse
import os
import sys



######################
## Get user's input ##
######################

## 1. Define parser ##
### Mandatory arguments
parser = argparse.ArgumentParser(description='Analyse fold-back inversion artefacts')
parser.add_argument('bam', help='Input bam file')
parser.add_argument('chromosome', help='Chromosome name. E.g. chr12')
parser.add_argument('sample_type', help='TUMOUR or NORMAL')
parser.add_argument('sample_id', help='Sample name')
parser.add_argument('output_dir', help='Output directory')


## 2. Parse user's input and initialize variables ##
args = parser.parse_args()
bam = args.bam
ref = args.chromosome
sample = args.sample_type
sample_id = args.sample_id
output_dir = args.output_dir
scriptName = os.path.basename(sys.argv[0])
version='1.0'

print()
print('***** ', scriptName, version, 'configuration *****')
print('*** Arguments ***')
print('bam file: ', bam)
print('chromosome: ', ref)
print('sample type: ', sample)
print('sample id: ', sample_id)
print('output directory: ', output_dir, "\n")



## FUNCTIONS ##
def get_ref_lengths(bam):
    '''
    Make dictionary containing the length for each reference

	Input:
		1. bam: indexed BAM file
	
	Output:
		1. lengths: Dictionary containing reference ids as keys and as values the length for each reference
    '''
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')
    
    ## Make dictionary with the length for each reference (move this code into a function)
    refLengths = dict(list(zip(bamFile.references, bamFile.lengths)))
    
    ## Close bam file
    bamFile.close()
    
    return refLengths


def overlap(begA, endA, begB, endB):
    '''
    Check if two ranges overlap    

    Input:
        1. begA: interval A begin position
        2. endA: interval A end position
        3. begB: interval B begin position
        4. endB: interval B end position

    Output:
        1. boolean: True (overlap) and False (no overlap)
        2. overlapLen: number of overlapping base pairs
    '''    
    maxBeg = max([begA, begB])
    minEnd = min([endA, endB])

    # a) Overlap
    if maxBeg <= minEnd:
        boolean = True
        overlapLen = minEnd - maxBeg + 1

    # b) No overlap
    else:
        boolean = False
        overlapLen = 0

    return boolean, overlapLen


def alignment_query_coord(CIGAR, supp_strand):
    '''
    Compute start and end of the alignment from CIGAR string

    Input:
        1. CIGAR: CIGAR string

    Output:
        1. queryBeg: alignment start on query
        2. queryEnd: alignment end on query
    '''
    ## 1. Read CIGAR string using proper module
    cigarTuples = Cigar(CIGAR)
    cigarList = list(cigarTuples.items()) if supp_strand == "+" else list(reversed(list(cigarTuples.items())))

    # 2. Infer where the alignment starts based on cigar 
    start = cigarList[0]
    startOperation = start[1]
    queryBeg = int(start[0]) if startOperation in ('S') else 0

    ## 3. Iterate over the operations and compute the alignment length
    alignmentLen = 0

    for cigarTuple in cigarList:

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ### Update reference alignment length
        ## a) Operations consuming query
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op I, tag 1, insertion to the reference
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if operation in ('M', 'I', '=', 'X'):
            alignmentLen += length

    return queryBeg, queryBeg+alignmentLen


def alignment_ref_length(CIGAR):
    '''
    Compute alignment length from CIGAR string

    Input:
        1. CIGAR: CIGAR string

    Output:
        1. alignmentLen: alignment length on reference
    '''
    ## 1. Read CIGAR string using proper module
    cigarTuples = Cigar(CIGAR)

    ## 2. Iterate over the operations and compute the alignment length
    alignmentLen = 0
    
    for cigarTuple in list(cigarTuples.items()):

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ### Update reference alignment length
        ## a) Operations consuming reference
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op D, tag 2, deletion from the reference
        # - Op N, tag 3, skipped region from the reference
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if operation in ('M', 'D', 'N', '=', 'X'):
            alignmentLen += length

    return alignmentLen


def collect_reads(ref, binBeg, binEnd, bam, confDict):
    '''
    Collect reads supporting fold-back inversion events in a genomic bin from a bam file

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * minMAPQ        -> minimum mapping quality

    Output:
        1. invList: List of fold-back inversion events 
        2. nb_parsedReads: Total number of reads
    
    '''
    ## Initialize inv list and counter
    invList = []
    nb_parsedReads = 0

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### 1. Filter alignments based on different criteria:
        ## Unmapped reads   
        if alignmentObj.is_unmapped == True:
            continue

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality
        if (MAPQ < confDict['minMAPQ']):
            continue
    
        ## Filter out secondary and supplementary alignments 
        if alignmentObj.is_secondary or alignmentObj.is_supplementary:
            continue
        
        ## counter of total reads
        nb_parsedReads += 1
            
        ## 2. Collect INVs
        fb_inv = is_fb_inv(alignmentObj, confDict['minMAPQ'])

        if fb_inv:
            invList.append(fb_inv)

    ## Close 
    bamFile.close()

    # return sv candidates
    return invList, nb_parsedReads


def is_fb_inv(alignmentObj, minMAPQ):
    '''
    For a read alignment, check if the read supports a fold-back inversion pattern

    Input:
        1. alignmentObj: pysam read alignment object
        2. minMAPQ: minimum MAPQ threshold

    Output:
        1. fb_inv: FB_INV object (None if no fold-back inversion found)
    '''
    fb_inv = None

    ## if read has supplementary alignments
    if alignmentObj.has_tag('SA'):
        sa_tags = alignmentObj.get_tag('SA').split(";") # Pattern ['chr12,82187865,-,10S1601M203D1812S,60,485', '']
        
        # A. if only one supplementary alignment --> Easiest case
        if len(' '.join(sa_tags).split()) == 1:
            supplAlignment = ' '.join(sa_tags).split()[0]
            supp_chrom = supplAlignment.split(",")[0]

            # B. if supp in the same chrom as primary
            if supp_chrom == alignmentObj.reference_name:
                primary_strand = '-' if alignmentObj.is_reverse else '+'
                supp_strand = supplAlignment.split(",")[2]
                
                # C. if inverted alignments
                if len(list(set([primary_strand, supp_strand]))) == 2:

                    supp_cigar = supplAlignment.split(",")[3]
                    supp_queryBeg, supp_queryEnd = alignment_query_coord(supp_cigar, supp_strand)
                    supp_refBeg = int(supplAlignment.split(",")[1])
                    supp_refEnd = supp_refBeg + alignment_ref_length(supp_cigar)
                    primary_queryBeg, primary_queryEnd = alignment_query_coord(alignmentObj.cigarstring, primary_strand)

                    # D. if supp and primary overlap in reference
                    if overlap(alignmentObj.reference_start, alignmentObj.reference_end, supp_refBeg, supp_refEnd)[0]:

                        # E. if mapq supp aln > minMAPQ:
                        if int(supplAlignment.split(",")[5]) >= minMAPQ:
                            
                            invDict = {}
                            invDict['primary_queryBeg'] = alignmentObj.query_alignment_start
                            invDict['supp_queryBeg'] = supp_queryBeg
                            first = min(invDict, key=invDict.get)

                            if first == 'primary_queryBeg':
                                fb_inv = FB_INV(alignmentObj.query_name, 
                                            alignmentObj.reference_name, 
                                            alignmentObj.reference_start, 
                                            alignmentObj.reference_end, 
                                            primary_strand, 
                                            primary_queryBeg, 
                                            primary_queryEnd,
                                            supp_refBeg, 
                                            supp_refEnd,
                                            supp_strand,
                                            supp_queryBeg, 
                                            supp_queryEnd, 
                                            alignmentObj.query_length,
                                            "primaryFirst")
                                                        
                            elif first == 'supp_queryBeg':
                                fb_inv = FB_INV(alignmentObj.query_name, 
                                            alignmentObj.reference_name, 
                                            supp_refBeg, 
                                            supp_refEnd,
                                            supp_strand,
                                            supp_queryBeg, 
                                            supp_queryEnd,
                                            alignmentObj.reference_start, 
                                            alignmentObj.reference_end, 
                                            primary_strand, 
                                            primary_queryBeg, 
                                            primary_queryEnd,
                                            alignmentObj.query_length,
                                            "suppFirst")
                                                       
    return fb_inv



## CLASSES ##
class FB_INV():
    '''
    Fold-back inversion class
    '''

    def __init__(self, read_name, ref, alnA_start, alnA_end, alnA_strand, alnA_query_start, alnA_query_end, alnB_start, alnB_end, alnB_strand, alnB_query_start, alnB_query_end, read_len, first):
        '''
        '''
        self.read_name = read_name
        self.read_len = read_len
        self.ref = str(ref)
        self.first = str(first)

        ## First alignment attributes
        self.alnA_start = int(alnA_start)
        self.alnA_end = int(alnA_end)
        self.alnA_strand = str(alnA_strand)
        self.alnA_query_start = int(alnA_query_start)
        self.alnA_query_end = int(alnA_query_end)
        
        ## First alignment attributes
        self.alnB_start = int(alnB_start)
        self.alnB_end = int(alnB_end)
        self.alnB_strand = str(alnB_strand)
        self.alnB_query_start = int(alnB_query_start)
        self.alnB_query_end = int(alnB_query_end)

        ## Estimate distances bewteen bkps
        a = self.alnA_start if self.alnA_strand == "+" else self.alnA_end
        b = self.alnA_end if self.alnA_strand == "+" else self.alnA_start
        c = self.alnB_start if self.alnB_strand == "+" else self.alnB_end
        d = self.alnB_end if self.alnB_strand == "+" else self.alnB_start
        
        self.a_b = abs(a-b)
        self.c_d = abs(c-d)
        self.b_c = abs(b-c)
        self.a_d = abs(a-d)



#################
## CALL METHOD ##
#################

## A. get reference length
ref_lengths = get_ref_lengths(bam)
beg = 0
end = ref_lengths[ref]


## B. Get FB_INV events
confDict = {}
confDict['minMAPQ'] = 20
invList, nb_parsedReads = collect_reads(ref, beg, end, bam, confDict)


## C. Write output
## C.1. Write artifact attributes to output file
outFile = output_dir + '/' + sample_id + '_' + sample + '_' + ref + '_' + 'INV.tsv'
with open(outFile, 'w') as f:

    # write header
    header = ['#read_name', 'ref', 'start', 'end', 'strand', 'dist_ab', 'dist_cd', 'dist_bc', 'dist_ad']
    f.write('\t'.join(header) + "\n")

    # write artifact attributes to a file
    for invObj in invList:

        invObj_list = [str(i) for i in [invObj.read_name, invObj.ref, 
                                        invObj.alnA_start, invObj.alnA_end, invObj.alnA_strand, 
                                        invObj.a_b, invObj.c_d, 
                                        invObj.b_c, invObj.a_d]]
        f.write('\t'.join(invObj_list) + "\n")


## C.2. Write counts to summary file
outFile2 = output_dir + '/' + sample_id + '_' + sample + '_' + ref + '_' + 'INV.summary.tsv'
with open(outFile2, 'w') as f:

    # write header
    header = ['#donor_id', 'sample_id', 'ref', 'nb_parsedReads', 'nb_invs', 'nb_invArtifacts']
    f.write( '\t'.join(header) + "\n" )

    # estimate total number of reads
    nb_invs = len(invList)

    # estimate total number of inversion-like artifacts
    nb_invArtifacts = len([i for i in invList if i.a_d < 150])

    # write row
    row = [str(item) for item in [sample_id, sample, ref, nb_parsedReads, nb_invs, nb_invArtifacts]]
    f.write( '\t'.join(row) + "\n" )


