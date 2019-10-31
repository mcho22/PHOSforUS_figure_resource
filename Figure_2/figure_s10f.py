#!/usr/bin/env python

import matplotlib.pyplot as plot
import numpy

from matplotlib import rcParams
import matplotlib.font_manager as fm
import matplotlib.pyplot as plot
from scipy import stats

########

file_path = "/Users/cmdb/Desktop/phospho_analysis/191007_figurebase" # Change as appropriate

########

tar_file = open("pre_calculated_data/log_190826.txt")

tar_list = []

for line in tar_file:
    if len(line.rstrip('\n')) > 1:
        if line[0] == "[":
            ssl = line.rstrip('\n').split()
            tar_list.append([float(ssl[0][1:-1]), float(ssl[1]), float(ssl[2])])

#print len(tar_list), len(tar_list)/6
tar_six = [[], [], [], [], [], []]
tar_len = 835

for i_1 in range(len(tar_list)):
    tdiv = (i_1%6)
    tar_six[tdiv].append(tar_list[i_1])

#print len(tar_six[0]), len(tar_six[1]), len(tar_six[2]), len(tar_six[3]), len(tar_six[4]), len(tar_six[5])

tar_target = [[], []]
for i_2 in range(tar_len):
    tar_target[0].append(tar_six[3][i_2][2])
    tar_target[1].append(tar_six[5][i_2][2])

####

prop = fm.FontProperties(fname=file_path+'/External_Scripts/calibri/Calibri.ttf')
rcParams.update({'font.size': 15})

####

abin = stats.gaussian_kde(tar_target[0])
bbin = stats.gaussian_kde(tar_target[1])
cdiff = numpy.minimum(abin(numpy.arange(-0.001, 0.00101, 0.00001)), bbin(numpy.arange(-0.001, 0.00101, 0.00001)))

plot.figure(figsize = [8, 6])
plot.hist(tar_target[0], bins = numpy.arange(0, 1.01, 0.005), color = "#800000", alpha = 0.7)
plot.hist(tar_target[1], bins = numpy.arange(0, 1.01, 0.005), color = "#000080", alpha = 0.7)
plot.xlim([0, 1])
plot.ylim([0, 60])

plot.xlabel("COREX/eSCAPE dG, R-square score distribution", fontproperties=prop, fontsize = 15)
plot.ylabel("Probability density", fontproperties=prop, fontsize = 15)
plot.xticks(fontproperties=prop, fontsize = 15)
plot.yticks(fontproperties=prop, fontsize = 15)
plot.savefig("results/reproduced_figure_s10f.png")


