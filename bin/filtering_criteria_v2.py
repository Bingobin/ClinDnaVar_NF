#!/usr/bin/env python3

import sys
import csv
import numpy as np
from scipy.stats import fisher_exact

def read_maf(maf_file, maf_dict):
    with open(maf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Hugo_Symbol"):
                continue
            tmp = line.split("\t")
            pos = "\t".join([tmp[4], tmp[5], tmp[10], tmp[12]])
            if pos in maf_dict:
                maf_dict[pos] += 1
            else:
                maf_dict[pos] = 1

def dict_pon(maf_file, maf_dict):
    with open(maf_file, 'r') as f:
        maf_reader = csv.DictReader(f, delimiter='\t')
        for line in maf_reader:
            pos = f"{line['Chromosome']}\t{line['Start_Position']}\t{line['Reference_Allele']}\t{line['Tumor_Seq_Allele2']}";
            maf_dict[pos] = {}
            maf_dict[pos]['NS'] = line['NS']
            maf_dict[pos]['NS2'] = line['NS2']
            maf_dict[pos]['ALT_R'] = line['ALT_R']
            maf_dict[pos]['REF_R'] = line['REF_R']

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <ch_mut> <pon.collapse> <whitelist>")
        sys.exit(1)

    ch_mut = sys.argv[1]
    pon = sys.argv[2]
    whitelist = sys.argv[3]

    PON = {}
    WL = {}

    read_maf(whitelist, WL)
    dict_pon(pon, PON)

    with open(ch_mut, 'r') as f:
        header = f.readline().strip()
        print(f"{header}\tPON\tWhiteList\tFisher_P\tFILTER_C")
#        maf_reader = csv.DictReader(f, delimiter='\t')
#        header = "\t".join(maf_reader)
#        print(f"{header}\tPON\tWhiteList\tFILTER_C")
        for line in f:
            line = line.strip()
            tmp = line.split("\t")
#            pos = f"{line['Chromosome']}\t{line['Start_Position']}\t{line['Reference_Allele']}\t{line['Tumor_Seq_Allele2']}";
            pos = "\t".join([tmp[4], tmp[5], tmp[10], tmp[12]])

            ch_pon = PON.get(pos, {}).get('NS', 0)
            ch_wl = WL.get(pos, 0)
            ch_filter = []
            ns2 = PON.get(pos, {}).get('NS2', 0)
            fisher_p = 0
            if int(ns2) >= 2:
                ch_filter.append('NS2')
            if np.mean(np.array(tmp[41].split("|")).astype(int)) <  500:
                ch_filter.append('LDEP')
            if np.mean(np.array(tmp[48].split("|")).astype(int)) <  5:
                ch_filter.append('LALT')
            if (len(tmp[10]) > 20 or len(tmp[12]) > 20) and int(tmp[51])  == 1:
                ch_filter.append('MAID')
            if int(tmp[51])  < 4 and ch_wl == 0:
                ch_filter.append('LCAL')
            if float(tmp[50])  > 0.35:
                 ch_filter.append('HVAF')
            if int(ch_pon) > 0:
                alt_mut = np.mean(np.array(tmp[48].split("|")).astype(int));
                dep_mut = np.mean(np.array(tmp[41].split("|")).astype(int));
                alt_pon = int(PON[pos]['ALT_R']);
                dep_pon = int(PON[pos]['ALT_R']) + int(PON[pos]['REF_R']);
                contingency_table = np.array([[alt_mut, alt_pon], [dep_mut, dep_pon]])
                odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
                fisher_p = p_value
                adjp = p_value * 305839
                if adjp > 0.05:
                    ch_filter.append('ARTI')
            if len(ch_filter)  == 0:
                filter_c = "PASS"
            else:
                filter_c = ";".join(ch_filter)
            print(f"{line}\t{ch_pon}\t{ch_wl}\t{fisher_p}\t{filter_c}")
