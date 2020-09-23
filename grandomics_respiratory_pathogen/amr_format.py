#!/usr/bin/env python
# coding=utf-8
import sys,os,re
import argparse,csv

parser = argparse.ArgumentParser(description='nanopore reads analysis amr pipline')
parser.add_argument('--input_file', metavar='input_file', help='')
parser.add_argument('--data_type', dest="data_type", type=str.lower, choices=['contig','read'], required=True,help='')
opt = parser.parse_args()

dict_drug={}
fo = open("out.amr.tsv", "w")

if opt.data_type == "read":##基于比对的结果
	for line in open(opt.input_file):
		arr=line.strip().split("\t")
		if arr[0] == "ARO Term":
			continue
		
		confidence_reads = float(arr[9])##completely mapped reads
		coverage_percent = float(arr[12])##average percent coverage
		mapq = float(arr[14])##average mapq (completely mapped reads)
		
		confidence = ""
		if mapq >= 50 and confidence_reads >= 30 and coverage_percent >= 80:
			confidence = "high"
		elif mapq <= 30 and confidence_reads < 10 and coverage_percent <= 50:
			confidence = "low"
		else:
			confidence = "middle"
		
		ARO_Term = arr[0]
		species = arr[8]
		drugs = arr[24]
		
		for drug in drugs.strip().split(";"):
			drug = drug.strip()
			if drug not in dict_drug:
				dict_drug[drug] = {"gene":{}, "species_or_Mechanism":{}, "confidence":{}}
			dict_drug[drug]["gene"][ARO_Term] = ""
			dict_drug[drug]["species_or_Mechanism"][species] = ""
			dict_drug[drug]["confidence"][confidence] = ""
			
	fo.write("Drugs\tConfidence\tSpecies\tARO Term\n")
	writer = csv.writer(fo, delimiter='\t', dialect='excel')

if opt.data_type == "contig":##基于组装的结果
	for line in open(opt.input_file):
		arr=line.strip().split("\t")
		if arr[0] == "ORF_ID":
			continue
		Cut_Off = arr[5]##Cut_Off
		if Cut_Off == "Perfect":
			confidence = "high"
		elif Cut_Off == "Strict":
			confidence = "middle"
		else:
			confidence = "low" ## Cut_Off==Loose
		
		ARO = arr[8]##Best_Hit_ARO
		drugs = arr[14]
		Resistance_Mechanism = arr[15]
		
		for drug in drugs.strip().split(";"):
			drug = drug.strip()
			if drug not in dict_drug:
				dict_drug[drug] = {"gene":{}, "species_or_Mechanism":{}, "confidence":{}}
			dict_drug[drug]["gene"][ARO] = ""
			dict_drug[drug]["species_or_Mechanism"][Resistance_Mechanism] = ""
			dict_drug[drug]["confidence"][confidence] = ""
		
	fo.write("Drugs\tConfidence\tResistance_Mechanism\tARO Term\n")
	writer = csv.writer(fo, delimiter='\t', dialect='excel')
		
for DRUG in dict_drug:
	confidence = ""
	confidence_list = dict_drug[DRUG]["confidence"].keys()
	if "high" in confidence_list:
		confidence = 'high'
	elif "middle" in confidence_list:
		confidence = "middle"
	elif "low" in confidence_list:
		confidence = "low"
	else:
		logger.info("confidence error...")
	writer.writerow([
					DRUG,
					confidence,
					(";").join(dict_drug[DRUG]["species_or_Mechanism"].keys()),
					(";").join(dict_drug[DRUG]["gene"].keys())
					])
fo.close()