#!/usr/bin/env python3
# Author: Zhang Siwen, zhangsiwen@grandomics.com
# History:
#     20200826, first version
#     20200921, change standard for positive decision
#     20201118, change positive/suspected/negative standard

import subprocess
import argparse
import pysam


# arguments passing
def GetArgs():
    parser = argparse.ArgumentParser(description='Filter bam, summary species info, reads and ratio.')
    require_args = parser.add_argument_group('required arguments')
    require_args.add_argument('--bam', dest='bam', help="bam file", required=True)
    require_args.add_argument('--bed', dest='bed', help="bed file with species and type info, chr\tstart\tend\tspecies\tname\ttype", required=True)
    require_args.add_argument('--out', dest='out', help='output prefix', required=True) 
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--depth', dest='depth', default=100, help="depth threshold, default 100", type=int)
    optional_args.add_argument('--cov', dest='cov', default=0.5, help="read-target overlap ratio, default 0.5", type=float)
    optional_args.add_argument('--mapq', dest='mapq', default=50, help="minimum mapping quality of reads to be kept,default 50", type=int)    
    args = parser.parse_args()
    return args


# record target region, subspecies relations
def ReadBedRelation(bed):
    relation = {} # species relation
    result = {} # species summary
    # target bed
    bedfile = open(bed)
    target = {}
    # store each region of bed
    for line in bedfile:
        try:
            ele = line.strip().split("\t")
            ref = ele[0]
            # check ref appear once or multiple times
            # e.g. target = {'NC_007373.1': [['4', '431']], 'NC_012564.1': [['2322', '2660'], ['4340', '5147']]}
            if ref in target:
                target[ref].append([ele[1],ele[2]])
            else:
                target[ref] = [[ele[1],ele[2]]]
            # ref name relation
            species = ele[3]
            # record relation of ref and species
            # e.g. relation = {'AC_000007.1': {'species': 'Adenovirus', 'type': 'Type2'}}
            relation[ref] = species
            relation[ref] = {'species': species, 'type':ele[5]}
            # set model of result
            result[species] = {'total_reads':0, 'total_ratio':0, 'total_cov':0, 'ref':[], 'type':[], 'reads':[], 'reads_ratio':[], 'target':[], 'coverage':[], 'cov_ratio':[]}
            # also add unclassified read in result
            #result['Unclassified'] = {'total_reads':0, 'total_ratio':0, 'total_cov':0, 'ref':['Unclassified'], 'type':['Unclassified'], 'reads':[], 'reads_ratio':[1], 'target':[0], 'coverage':[0], 'cov_ratio':[0]}
        except:
            continue
    bedfile.close()
    return target, relation, result


# filter bam file
def FilterBam(inbam, outbam, target, cov=0.5, mapq=50):
    samfile = pysam.AlignmentFile(inbam, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template = samfile)
    for record in samfile:
        # keep unmapped reads
        if record.is_unmapped:
            #outfile.write(record)
            continue
        # ignore qcfail, duplicate
        if record.is_qcfail or record.is_duplicate:
            continue
        # supplementary alignment
        try:
            SA=record.get_tag('SA')
        except:
            SA=''
        # ignore secondary or supplementary
        if record.is_secondary or SA:
            continue
        # ignore mapping_quality fail
        if record.mapping_quality < mapq:
            continue
        # calculate overlap ratio and write record
        for pos in target[record.reference_name]:
            start, end = int(pos[0]), int(pos[1])
            overlap = record.get_overlap(start, end) # start, end
            # keep record with overlap ratio > threshold
            if overlap and overlap/record.query_length >= cov and overlap/(end-start) >= cov:
                outfile.write(record)
    samfile.close()
    outfile.close()
    subprocess.run(["samtools", "index", outbam])


# calculate reads and coverage for each ref
def BamReadsCoverage(outbam, target, relation, result, depth=100):
    outfile = pysam.AlignmentFile(outbam, "rb")
    # calculate read number and coverage for each subspecies and merge it to corresponding species
    for ref in outfile.references:
        # ignore non-bed regions
        if ref not in target:
            continue
        # read number mapped to ref
        read = 0
        # total length of ref target
        length = 0
        coverage = []
        # check each region in bed
        for pos in target[ref]:
            start, end = int(pos[0]), int(pos[1])
            # read number
            read += outfile.count(ref, start, end)
            length += end-start
            # four base coverage
            basecov = outfile.count_coverage(ref, start, end)
            # total depth for each pos
            totalcov = [basecov[0][i]+basecov[1][i]+basecov[2][i]+basecov[3][i] for i in range(len(basecov[0]))]
            # [list] count pos with depth higher than threshold
            coverage += [k for k in totalcov if k > depth] 
        
        # add species info to result
        # species name of ref
        species = relation[ref]['species']
        # ref name, e.g. NC_002207.1
        result[species]['ref'].append(ref)
        # subtype name, e.g. Lee-seg4
        result[species]['type'].append(relation[ref]['type'])
        # reads number of subspecies, list of int()
        result[species]['reads'].append(read)
        # length of subspecies, list of int()
        result[species]['target'].append(length)
        # base number of subspecies with depth > threshold, list of int()
        result[species]['coverage'].append(len(coverage))
        # coverage ratio of subspecies, list of float(2)
        result[species]['cov_ratio'].append(round(len(coverage)/length*100,2))
    
    # also record unmapped reads
    #result['Unclassified']['reads'].append(outfile.unmapped)
    outfile.close()
    
    # total reads
    total_reads = sum([i for k,v in result.items() for i in v['reads']])
    # get total info for each species
    for ref in result:
        result[ref]['total_reads'] = sum(result[ref]['reads'])
        result[ref]['total_ratio'] = round(result[ref]['total_reads']/total_reads*100,2)
        # choose coverage of max read number as total coverage
        idx = result[ref]['reads'].index(max(result[ref]['reads']))
        result[ref]['total_cov'] = result[ref]['cov_ratio'][idx]
        # calculate reads ratio for each subspecies
        if result[ref]['total_reads'] != 0:
            result[ref]['reads_ratio'] = [round(r/result[ref]['total_reads']*100,2) for r in result[ref]['reads']]
        else:
            # if total_reads == 0, set [0] with length of total_reads
            result[ref]['reads_ratio'] = [0]*len(result[ref]['reads'])
    return result, total_reads


# output summary file
def WriteOutSummary(outfile, result, total_reads, depth):
    # store subspecies infomation
    info = {}
    out = open(outfile,"w")
    # write header
    out.write("\t".join(['Species', 'Total reads', 'Reads ratio', str(depth)+'X coverage%', 'Positive', 'info(name,reads,ratio,coverage)'])+"\n")
    for species in result:
        info[species] = []
        # join info of subspecies for each species
        for i in range(len(result[species]['ref'])):
            # e.g. Lee-seg4(NC_002207.1,0,0%,0.0%);Lee-seg6(NC_002209.1,0,0%,0.0%);Lee-seg7(NC_002210.1,0,0%,0.0%)
            info[species].append(result[species]['type'][i]+"("+result[species]['ref'][i]+","+str(result[species]['reads'][i])+","+str(result[species]['reads_ratio'][i])+"%,"+str(result[species]['cov_ratio'][i])+"%)")
        # positive decision, reads_ratio > 10 and cov_depth > 90
        #if result[species] == "Unclassified": # modified 20200921
            #decide = "No"
        #el
        if (result[species]['total_ratio'] >= 10 or result[species]['total_reads'] >= 1000) and (result[species]['total_cov'] >= 30): # modified 20201118
            decide = "Yes"
        elif (result[species]['total_ratio'] >= 1 or result[species]['total_reads'] >= 100) and (result[species]['total_cov'] >= 10): # modified 20201118
            decide = "Suspected"
        elif result[species]['total_cov'] <= 10: # modified 20201118
            decide = "No" # modified 20201118
        elif result[species]['total_ratio'] <= 1 and result[species]['total_reads'] <= 100: # modified 20201118
            decide = "No" # modified 20201118
        else:
            decide = "Warning" # modified 20200921
        # write out species stats
        out.write(species+"\t"+str(result[species]['total_reads'])+"\t"+str(result[species]['total_ratio'])+"%\t"+str(result[species]['total_cov'])+"%\t"+decide+"\t")
        # subspecies info
        out.write(";".join(info[species])+"\n")
        # e.g. 
        # Species	Total reads	Reads ratio	100X coverage%	Positive	info(name,reads,ratio,coverage)
        # Influenza B virus	0	0.0%	0.0%	no	Lee-seg4(NC_002207.1,0,0%,0.0%);Lee-seg6(NC_002209.1,0,0%,0.0%);Lee-seg7(NC_002210.1,0,0%,0.0%)
    # finally write out Total reads info
    out.write("Total reads\t"+str(total_reads)+"\t100%\t-\t-\t-\n")
    if total_reads < 1000:
        out.write("Attention: Reads of this sample is not enough, the cutoff is 1000, now you only have "+str(total_reads)+".\n")
        out.write("Attention: The positive results are not warranted, re-sequencing of this sample is recommanded.\n")
    out.close()


def main():
    args = GetArgs()
    # record target region, subspecies relations
    target, relation, result = ReadBedRelation(args.bed)

    # filter bam file
    FilterBam(args.bam, args.out+".bam", target, args.cov, args.mapq)

    # calculate reads and coverage for each ref
    result, total_reads = BamReadsCoverage(args.out+".bam", target, relation, result, args.depth)
    
    # write out summary
    WriteOutSummary(args.out+".summary.tsv", result, total_reads, args.depth)


if __name__=="__main__":
    main()
    
