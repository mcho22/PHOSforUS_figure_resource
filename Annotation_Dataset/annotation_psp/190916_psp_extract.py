#!/usr/bin/env python

import datetime
import numpy
import sys

####

label_site = []
label_seq = []
sitelist = []
seqlist = []

det_1 = 0
det_2 = 0

with open("Phosphorylation_site_dataset") as file_1:
    for line in file_1:
        ssl = line.rstrip('\n').split('\t')
        if det_1 > 0:
            sitelist.append(ssl) # 2, 4, 6, 9, 10, 11, 12, 13
        if len(ssl) > 2:
            if ssl[0] == "GENE":
                label_site.extend(ssl)
                det_1 += 1

with open("Phosphosite_PTM_seq.fasta") as file_2:
    for line in file_2:
        ssl = line.rstrip('\n').split('|')
        if det_2 >= 3:
            if line[0] == ">":
                label_seq.append([ssl[2], ssl[3]])
                seqlist.append("")
            else:
                seqlist[-1] += line.rstrip('\n')
        det_2 += 1

####

human_sites = []
human_seq = []

for i_1 in range(len(sitelist)):
    if sitelist[i_1][6] == "human":
        human_sites.append(sitelist[i_1])

for i_2 in range(len(label_seq)):
    if label_seq[i_2][0] == "human":
        human_seq.append([label_seq[i_2][0], label_seq[i_2][1], "______________" + seqlist[i_2] + "______________"])

# Total 371321 sites // Human 239668 sites
# Total 57052 sequences // Human 25652 sequences

####

extend_sites = [[], [], []]

for i_3 in range(len(human_sites)):
    tar_site = int(human_sites[i_3][4][1:-2]) + 13
    tar_type = -1
    if len(human_sites[i_3][10]) > 0:
        tar_type += 1
    elif len(human_sites[i_3][11]) > 0:
        tar_type += 2
    else:
        tar_type += 3
    site_type = -1
    if human_sites[i_3][9][7] == "s":
        site_type += 1
        if human_sites[i_3][9][8] == "P":
            site_type += 3
    elif human_sites[i_3][9][7] == "t":
        site_type += 2
        if human_sites[i_3][9][8] == "P":
            site_type += 3
    elif human_sites[i_3][9][7] == "y":
        site_type += 3
    for i_4 in range(len(human_seq)):
        if human_sites[i_3][2] == human_seq[i_4][1]:
            if human_sites[i_3][9] == human_seq[i_4][2][tar_site-7:tar_site+8]:
                extend_sites[tar_type].append([human_sites[i_3][2], str(tar_site), str(site_type), human_seq[i_4][2][tar_site-14:tar_site+15].upper()])
                break
    if i_3 % 10000 == 9999:
        print i_3, len(extend_sites[0]), len(extend_sites[1]), len(extend_sites[2]), len(extend_sites[0]) + len(extend_sites[1]) + len(extend_sites[2]), str(datetime.datetime.now())[8:19]

# Human-LTP-LT 13520, Human-MS-LT 139887, Human-MS-CST 86261 (sorted)
# Human-LTP-LT 13492, Human-MS-LT 139733, Human-MS-CST 86175 (after-re-matching, 268 sites lost)

print len(human_sites), len(extend_sites[0]), len(extend_sites[1]), len(extend_sites[2]), len(extend_sites[0]) + len(extend_sites[1]) + len(extend_sites[2]), len(human_sites)-len(extend_sites[0]) + len(extend_sites[1]) + len(extend_sites[2]), str(datetime.datetime.now())[8:19]

fname_list = ["psp_ltp_", "psp_ms_lit_", "psp_ms_cst_"]

for i_5 in range(3):
    for i_6 in range(5):
        out_file = open(fname_list[i_5] + str(i_6) + ".txt", 'w')
        out_count = 0
        for i_7 in range(len(extend_sites[i_5])):
            if int(extend_sites[i_5][i_7][2]) == i_6:
                sst = extend_sites[i_5][i_7][0] + '\t' + extend_sites[i_5][i_7][1] + '\t' + extend_sites[i_5][i_7][3] + '\n'
                out_file.write(sst)
                out_count += 1
        print i_5, i_6, out_count
        out_file.close()

