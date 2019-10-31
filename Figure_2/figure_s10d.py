#!/usr/bin/env python

import datetime
import gc
import matplotlib.pyplot as plot
import numpy
import os
import subprocess
import sys

from Bio.SubsMat import MatrixInfo
from scipy.stats import pearsonr
from sklearn import metrics
from sklearn.linear_model import LinearRegression

from matplotlib import rcParams
import matplotlib.font_manager as fm
import matplotlib.pyplot as plot

########

file_path = "/Users/cmdb/Desktop/phospho_analysis/191007_figurebase" # Change as appropriate

########

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def site_pairwise(aa_list, matrix, gap_s):
    sp_score = []
    for sp_1 in range(len(aa_list)):
        max_score = float(score_match((aa_list[sp_1], aa_list[sp_1]), matrix))
        for sp_2 in numpy.arange(sp_1+1, len(aa_list), 1):
            sp_pair = (aa_list[sp_1], aa_list[sp_2])
            if 'X' in sp_pair:
                sp_score.append(gap_s / max_score)
            else:
                sp_score.append(float(score_match(sp_pair, matrix)) / max_score)
    return numpy.mean(sp_score)

########

tar_file = open(file_path+"/Annotation_Dataset/oma_tf_alignment/0616_P04150.fasta")
tar_label = []
tar_list = []
human_list = []
det_x = 0

for line in tar_file:
    sst = line.rstrip('\n')
    if line[0] == ">":
        if line[0] == ">": # line[0] == ">" // line[1:6] in species_sel
            tar_label.append(sst[1:])
            tar_list.append("")
            det_x += 1
        else:
            det_x = 0
    else:
        if det_x > 0:
            for i_0 in range(len(sst)):
                if sst[i_0] == "-":
                    tar_list[-1] += "X"
                else:
                    tar_list[-1] += sst[i_0]

for h_1 in range(len(tar_label)):
    if tar_label[h_1][0:5] == "HUMAN":
        for h_2 in range(len(tar_list[h_1])):
            human_list.append(tar_list[h_1][h_2])

tar_file.close()

print len(tar_label), len(human_list)

########

xsites = []
sequence_list = []
value_list = [[], [], [], [], [], [], [], [], [], []]

for i_1 in range(len(tar_label)): #len(tar_label)
    ####
    temp_long = []
    temp_short = []
    escape_list = []
    ####
    temp_file = open("input_temp/"+tar_label[i_1]+".fasta", 'w')
    temp_file.write(">"+tar_label[i_1])
    temp_file.write('\n')
    temp_file.write(tar_list[i_1])
    temp_file.write('\n')
    temp_file.close()
    #### eSCAPE prediction
    subprocess.check_output(["perl", "./eScape_predictor_REP.pl", "input_temp", "output_escape"])
    escape_open = open("output_escape/"+tar_label[i_1]+".eScapev8_pred")
    with open("output_escape/"+tar_label[i_1]+".eScapev8_pred") as escape_open:
        for line in escape_open:
            escape_list.append([float(i) for i in line.rstrip('\n').split()])
    #### Temporary file clear
    os.remove("input_temp/"+tar_label[i_1]+".fasta")
    os.remove("output_escape/"+tar_label[i_1]+".eScapev8_pred")
    os.remove("output_escape/"+tar_label[i_1]+".eScapev8_log")
    gc.collect()
    #### Data output
    escape_transpose = numpy.transpose(escape_list)
    sequence_sub = []
    for i_2 in range(len(tar_list[i_1])):
        sequence_sub.append(tar_list[i_1][i_2])
        if i_1 == 0:
            xsites.append(i_2+1)
    sequence_list.append(sequence_sub)
    for i_2_1 in range(2):
        value_list[i_2_1].append(escape_transpose[i_2_1])
    for i_3 in range(8):
        value_list[2+i_3].append(escape_transpose[i_3])
    if i_1 % 10 == 9:
        print "-- Sequence #"+str(i_1+1)+" calculation completed: "+tar_label[i_1]

####
################################################################
####

x_sel = []
seq_sel = []
value_sel = [[], [], [], [], [], [], [], [], [], []]

for i_4 in range(len(sequence_list)):
    seq_sel.append([])
    for i_5 in range(len(value_list)):
        value_sel[i_5].append([])

for i_7 in range(len(xsites)):
    if human_list[i_7] != "X":
        x_sel.append(i_7+1)
        for i_8 in range(len(seq_sel)):
            seq_sel[i_8].append(sequence_list[i_8][i_7])
            for i_9 in range(10):
                value_sel[i_9][i_8].append(value_list[i_9][i_8][i_7])

#print len(x_sel), len(seq_sel[0]), len(value_sel[0][0])
seq_array = numpy.asarray(seq_sel)
value_array = numpy.asarray(value_sel)

####
################################################################
####

blosum = MatrixInfo.blosum62
sequence_dscore = []
for i_4 in range(len(x_sel)):
    sequence_dscore.append(site_pairwise(seq_array[:, i_4], blosum, -1.0))

sequence_ascore = []
windowbasis = 5
windowsize = windowbasis*2-1

for i_5 in range(len(sequence_dscore)-windowsize+1):
    sequence_ascore.append(numpy.mean(sequence_dscore[i_5:(i_5+windowsize)]))

print "## Sequence similarity calculation completed"

####
################################################################
####

####
################################################################
####

rt = 1.0
rt_rx = [1.0, 1.0, rt, rt, rt, rt, rt, rt, rt, rt]
value_average = [[], [], [], [], [], [], [], [], [], []]
value_processed = [[], [], [], [], [], [], [], [], [], []]

for i_8 in range(len(sequence_dscore)-windowsize+1):
    for i_9 in range(3):
        vsubmean = []
        for i_10 in range(len(seq_sel)):
            vsubmean.append(numpy.mean(value_array[i_9, i_10, i_8:(i_8+windowsize)])/rt_rx[i_9])
        vmean = numpy.mean(vsubmean)
        value_average[i_9].append(vmean)
        vstdv = numpy.std(vsubmean)
        #value_processed[i_9].append(1.0-((vstdv*rt_rx[i_9])/numpy.fabs(numpy.mean(value_array[i_9, :, :]))))
        #vstdv = numpy.mean([1.0-((numpy.fabs(vsubmean[i]-vmean)*rt_rx[i_9])/ numpy.fabs(numpy.mean(value_array[i_9, :, :]))) for i in range(len(vsubmean))])
        value_processed[i_9].append(vstdv)
    if i_8 % 100 == 99:
        print "-- BC#:", i_8+1

####
################################################################
####

sa_mean = numpy.mean(sequence_ascore)
sa_stdv = numpy.std(sequence_ascore)

dg_mean = numpy.mean(value_processed[2])
dg_stdv = numpy.std(value_processed[2])

sq_sd = [(float(sequence_ascore[i]-sa_mean)/sa_stdv) for i in range(len(sequence_ascore))]
dg_sd = [(float(dg_mean-value_processed[2][k])/dg_stdv) for k in range(len(sequence_ascore))]

r4 = LinearRegression()
r4.fit(numpy.asarray(value_average[2]).reshape((-1, 1)), sq_sd)
r4_fit = r4.predict(numpy.asarray(value_average[2]).reshape((-1, 1)))
#print r4.coef_, r4.intercept_, metrics.r2_score(sq_sd, r4_fit)

r6 = LinearRegression()
r6.fit(numpy.asarray(value_average[2]).reshape((-1, 1)), dg_sd)
r6_fit = r6.predict(numpy.asarray(value_average[2]).reshape((-1, 1)))
#print r6.coef_, r6.intercept_, metrics.r2_score(dg_sd, r6_fit)

plot.figure(figsize = [8, 6])
prop = fm.FontProperties(fname=file_path+'/External_Scripts/calibri/Calibri.ttf')
rcParams.update({'font.size': 15})

plot.plot(value_average[2], sq_sd, linestyle = '', marker = 'o', markersize = 15, color = "#B00000", alpha = 0.05)
plot.plot(value_average[2], dg_sd, linestyle = '', marker = 'o', markersize = 15, color = "#0000B0", alpha = 0.05)
plot.plot(value_average[2], r4_fit, linestyle = ':', color = '#B00000')
plot.plot(value_average[2], r6_fit, linestyle = ':', color = '#0000B0')

plot.xlabel("COREX/eSCAPE dG (kcal/mol)", fontproperties=prop, fontsize = 15)
plot.ylabel("Normalized conservation scores", fontproperties=prop, fontsize = 15)
plot.xticks(fontproperties=prop, fontsize = 15)
plot.yticks(fontproperties=prop, fontsize = 15)
plot.savefig("results/reproduced_figure_s10d.png")


