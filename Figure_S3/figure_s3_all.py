#!/usr/bin/env python

####

import datetime
import matplotlib.pyplot as plot
import numpy
import random
import sys

from scipy.stats import gaussian_kde

####

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_", "#", "$", "&"]

nhp_raw = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, 0.0, -0.8, -0.7, -1.3]
thp_raw = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, 0.0, -4.5, -4.5, -4.5]
nhp_list = [((i/9.0)+0.5) for i in nhp_raw]
thp_list = [((i/9.0)+0.5) for i in thp_raw]

neg_nch = [0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
neg_tch = [0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, -2.0]
pos_nch = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

calc_range = 10

####

def hp_wave(hp_crange):
    hw_len = hp_crange + hp_crange + 1
    hw_half = hp_crange
    hw_weight = []
    for hw_1 in range(hw_len):
        if hw_1 == 0:
            hw_weight.append(1.0)
        elif hw_1 > 0 and hw_1 <= hw_half:
            hw_weight.append(hw_1+1.0)
        elif hw_1 > hw_half:
            hw_weight.append(hw_weight[-1]-1.0)
    return hw_weight

hp_weight = hp_wave(calc_range)

####

vlist = [[[], [], [], [], []], [[], [], [], [], []]]

with open("phospho_structure_human.txt") as tar_file:
    for line in tar_file:
        ssl = line.rstrip('\n').split('\t')
        st_stat = int(ssl[3])
        ph_stat = 0
        if ssl[2][14] == "#":
            ph_stat += 0
            if ssl[2][15] == "P":
                ph_stat += 3
        elif ssl[2][14] == "$":
            ph_stat += 1
            if ssl[2][15] == "P":
                ph_stat += 3
        elif ssl[2][14] == "&":
            ph_stat += 2
        vlist[0][ph_stat].append(ssl[2])

for i_1 in numpy.arange(0,1,1):
    for i_2 in range(5):
        tar_values = [[], [], [], [], [], [], [], []]
        for i_3 in range(len(vlist[i_1][i_2])):
            numeric = []
            for i_4 in range(29):
                det_1 = 0
                for i_5 in range(24):
                    if vlist[i_1][i_2][i_3][i_4] == aa_list[i_5]:
                        numeric.append(i_5)
                        det_1 += 1
                if det_1 < 1:
                    numeric.append(20)
                    det_1 += 1
            neg_pre = []
            neg_post = []
            pos_list = []
            hp_pre = []
            hp_post = []
            for i_6 in range(29):
                neg_pre.append(neg_nch[numeric[i_6]])
                neg_post.append(neg_tch[numeric[i_6]])
                pos_list.append(pos_nch[numeric[i_6]])
                hp_pre.append(nhp_list[numeric[i_6]])
                hp_post.append(thp_list[numeric[i_6]])
            sub_v = [numpy.mean(neg_pre[(14-calc_range):(15+calc_range)]), numpy.mean(neg_post[(14-calc_range):(15+calc_range)]), numpy.mean(pos_list[(14-calc_range):(15+calc_range)]), numpy.vdot(hp_pre[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(hp_post[(14-calc_range):(15+calc_range)], hp_weight)/121.0]
            tar_values[0].append(sub_v[0])
            tar_values[1].append(sub_v[1])
            tar_values[2].append(sub_v[2])
            tar_values[3].append(sub_v[2]-sub_v[0])
            tar_values[4].append(sub_v[2]-sub_v[1])
            tar_values[5].append(sub_v[3])
            tar_values[6].append(sub_v[4])
        
        print numpy.mean(tar_values[3]), numpy.mean(tar_values[4]), numpy.mean(tar_values[5]), numpy.mean(tar_values[6])
        
        X_n, Y_n = numpy.mgrid[0:1:201j, 0:1:200j]
        positions = numpy.vstack([X_n.ravel(), Y_n.ravel()])
        v_1 = numpy.vstack([tar_values[3], tar_values[4]])
        k_1 = gaussian_kde(v_1)
        Z_1 = numpy.reshape(k_1(positions).T, X_n.shape)
        
        plot.figure(figsize = [7, 7])
        #plot.contour(numpy.rot90(numpy.fliplr(Z_1)), cmap = "Blues", levels = [1.0, 3.0, 9.0], extent = [0, 1, 0, 1])
        plot.imshow(numpy.rot90(Z_1), cmap = "Blues", extent = [0, 1, 0, 1])
        p_types = ["b", "d", "e", "a", "c"]
        
        plot.plot([0, 1], [0.25, 0.25], 'k--')
        plot.plot([0, 1], [0.35, 0.35], 'k--')
        plot.plot([0.25, 0.25], [0, 1], 'k--')
        plot.plot([0.35, 0.35], [0, 1], 'k--')
        plot.xlabel("Pre-phosphorylation")
        plot.ylabel("Post-phosphorylation")
        plot.legend(loc = 4)
        plot.savefig("reproduced_figure_s3"+p_types[i_2]+".png")
        plot.clf()


