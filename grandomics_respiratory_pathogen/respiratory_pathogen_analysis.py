#!/usr/bin/env python3
# Author: Zhang Siwen, zhangsiwen@grandomics.com
# History:
#     20200827, first version
#     20200903, add barcode-sample relation file 
#     20200921, change Demultiplex(), let user choose if the data need demultiplexing


import os
import argparse
import subprocess
import multiprocessing

# arguments passing
def GetArgs():
    parser = argparse.ArgumentParser(description='Respiratory Pathogen Analysis.')
    require_args = parser.add_argument_group('required arguments')
    require_args.add_argument('--cell', dest='cell', help="cell number", required=True)
    require_args.add_argument('--kit', dest='kit', help="barcode kits", required=True)
    require_args.add_argument('--ref', dest='ref', help="reference path", required=True)
    require_args.add_argument('--tool', dest='tool', help="tool script bam_filter_stats_respiratory_pathogen.py path", required=True)
    require_args.add_argument('--bed', dest='bed', help="bed file with species and type info, chr\tstart\tend\tspecies\tname\ttype", required=True)
    require_args.add_argument('--relation', dest='relation', help="csv file of barcode and sample relation, 'barcode01,sample1'", required=True)
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--threads', dest='threads', default='4', help="threads, default 4")
    optional_args.add_argument('--process', dest='process_num', default='4', help="number of sub-process, default 4")
    optional_args.add_argument('--depth', dest='depth', default='100', help="depth threshold, default 100")
    optional_args.add_argument('--cov', dest='cov', default='0.5', help="overlap ratio, default 0.5")
    optional_args.add_argument('--mapq', dest='mapq', default='50', help="minimum mapping quality of reads to be kept,default 50")
    
    args = parser.parse_args()
    return args



# sub-process control
def RunSubprocess(c):
    pid = subprocess.Popen(c, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while pid.poll() == None: #检查子进程是否已被终止。设置并返回 returncode 属性。否则返回 None。
        out = pid.stdout.readline().strip()
        if out:
            print("subprocess output: ", out)


# multi-process control
def RunMultiprocess(commands, process_num, threads):
    # calculate a reasonable parallel job number
    process_num = min(int(os.cpu_count()/int(process_num)), int(threads))
    
    # create nultiprocess pool
    multiprocess_pool = multiprocessing.Pool(processes=process_num)
    for c in commands:
        multiprocess_pool.apply_async(RunSubprocess, args=(c,)) 
    print('Waiting for all subprocesses done...')
    multiprocess_pool.close()
    multiprocess_pool.join()
    print('All subprocesses done.')


# parse sample barcode relation file
def SampleBarcodeMap(relation): # added 20200921
    # record barcode-sample name relation 
    samplemap = {}
    relationfile = open(relation)
    for line in relationfile:
        ele = line.split(",")
        #ele = line.split("\t")
        try:
            samplemap[ele[0].strip()] = ele[1].strip()
        except:
            continue
    relationfile.close()

    return samplemap


# demultiplex
def Demultiplex(cell, kit, threads, process_num):
    # get system env
    basecalled_dir = os.getenv('BASECALLED_DIR')
    
    # check if basecalled_dir exists
    if not basecalled_dir: 
        raise Exception("Can not find $BASECALLED_DIR in Environment variables")
    
    # assign cell dir name 
    cell_dir = os.path.join(basecalled_dir, cell, "fastq_pass")
    
    # check if cell dir exists
    if not os.path.isdir(cell_dir): 
        raise Exception("Can not find fastq_pass/ directory in cell %s in $BASECALLED_DIR" % cell_dir)
    
    # decide if demultiplex is needed # added 20200921
    barcodes = {}
    if kit != "none":
        # run guppy_barcoder to demultiplex
        c = "/opt/conda/envs/galaxy/ont-guppy-cpu/bin/guppy_barcoder -i %s -s demultiplex --barcode_kits \"%s\" --trim_barcodes -t %s -q 0 -x auto --num_extra_bases_trim 0" % (cell_dir, kit, threads) # add '-x auto --num_extra_bases_trim 10' 20200923
        RunSubprocess(c)
        
        # find output folder name of each barcode
        dirs =["demultiplex/"+i for i in [d for d in os.listdir("demultiplex") if os.path.isdir('/'.join(["demultiplex",d]))]]
        fastqs = ['/'.join([i,os.listdir(i)[0]]) for i in dirs]
        
        # record barcode and file relation
        for f in fastqs:
            bc = f.split("/")[1]
            barcodes[bc] = f
        
    else: # added 20200921
        # demultiplexing already done
        os.system("mkdir demultiplex")
        commands = ["cat "+cell_dir+"/"+bc+"/*fastq > demultiplex/"+bc+".fastq" for bc in os.listdir(cell_dir)]
        RunMultiprocess(commands, process_num, threads)
        
        # record barcode and file relation
        for bc in os.listdir(cell_dir):
            barcodes[bc] = "demultiplex/"+bc+".fastq"
        
    return barcodes
    

# mapping
def MappingAndSummary(barcodes, ref, threads, process_num, tool, bed, relation, depth, cov, mapq):
    # record barcode-sample name relation 
    samplemap = SampleBarcodeMap(relation)
    
    # give warning if barcode in relation file but fastq file not found # added 20200921
    for bc in samplemap:
        if bc not in barcodes.keys():
            print("Warning: " + bc + " for " + samplemap[bc] + " not found. Continue...")
        else:
            print("Continue: " + bc + " for " + samplemap[bc] + " found. Continue...")
    
    # call minimap2
    #os.system("mkdir output")
    commands = ["minimap2 --MD -L -Y -t %s --secondary=no -ax map-ont %s %s|samtools view -b - |samtools sort -m 4G -o %s.bam - && samtools index %s.bam && python3 %s --bam %s.bam --bed %s --out %s.filter --depth %s --cov %s --mapq %s" % (threads, ref, fq, samplemap[bc], samplemap[bc], tool, samplemap[bc], bed, samplemap[bc], depth, cov, mapq) for bc,fq in barcodes.items() if bc in samplemap.keys()]
    
    RunMultiprocess(commands, process_num, threads)
    
    
def main():
    args = GetArgs()
    barcodes = Demultiplex(args.cell, args.kit, args.threads, args.process_num)
    MappingAndSummary(barcodes, args.ref, args.threads, args.process_num, args.tool, args.bed, args.relation, args.depth, args.cov, args.mapq)
    print("Analysis finished.")
    

if __name__=="__main__":
    main()
    
