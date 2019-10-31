#!/usr/bin/env python

####

import datetime
import matplotlib.pyplot as plot
import numpy
import random
import sys

from scipy.stats import gaussian_kde, ttest_ind, f_oneway, kruskal

import matplotlib.cm as cm
import matplotlib.font_manager as fm
from matplotlib import rcParams

####

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_", "#", "$", "&"]

#nhp_raw = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, 0.0, -0.8, -0.7, -1.3]
#thp_raw = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, 0.0, -4.5, -4.5, -4.5]
#nhp_list = [((i/9.0)+0.5) for i in nhp_raw]
#thp_list = [((i/9.0)+0.5) for i in thp_raw]

nhp_list = [0.37, 0.25, 0.30, 0.42, 0.17, 0.13, 0.20, 0.39, 0.56, 0.24, 0.36, 0.27, 1.00, 0.53, 0.38, 0.24, 0.32, 0.39, 0.25, 0.25, 0.0, 0.24, 0.32, 0.25] ## PII-pre
thp_list = [0.37, 0.25, 0.30, 0.42, 0.17, 0.13, 0.20, 0.39, 0.56, 0.24, 0.36, 0.27, 1.00, 0.53, 0.38, 0.24, 0.32, 0.39, 0.25, 0.25, 0.0, 0.26, 0.90, 0.16] ## PII-phos

neg_nch = [0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
neg_tch = [0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, -2.0]
pos_nch = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#target_pair = [[1, 0], [4, 3], [3, 0], [4, 1]]

calc_range = int(sys.argv[1])

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
verr = 0

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
        #vlist[st_stat][ph_stat].append(ssl[2])
        sst = ""
        for i_0 in range(29):
            if i_0 != 14:
                if ssl[2][i_0] == "#":
                    sst += "S"
                elif ssl[2][i_0] == "$":
                    sst += "T"
                elif ssl[2][i_0] == "&":
                    sst += "Y"
                else:
                    sst += ssl[2][i_0]
            else:
                sst += ssl[2][i_0]
        vlist[st_stat][ph_stat].append(sst)

print len(vlist[0][0]), len(vlist[0][1]), len(vlist[0][2]), len(vlist[0][3]), len(vlist[0][4])
print len(vlist[1][0]), len(vlist[1][1]), len(vlist[1][2]), len(vlist[1][3]), len(vlist[1][4])

"""    
for i_1 in range(5):
    tar_values = [[[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []]]
    for i_2 in range(2):
        for i_3 in range(len(vlist[i_2][i_1])):
            numeric = []
            for i_4 in range(29):
                det_1 = 0
                for i_5 in range(24):
                    if vlist[i_2][i_1][i_3][i_4] == aa_list[i_5]:
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
            tar_values[i_2][0].append(sub_v[0])
            tar_values[i_2][1].append(sub_v[1])
            tar_values[i_2][2].append(sub_v[2])
            tar_values[i_2][3].append(sub_v[2]-sub_v[0])
            tar_values[i_2][4].append(sub_v[2]-sub_v[1])
            tar_values[i_2][5].append(sub_v[3])
            tar_values[i_2][6].append(sub_v[4])
        
    print numpy.mean(tar_values[0][3]), numpy.mean(tar_values[0][4]), numpy.mean(tar_values[0][5]), numpy.mean(tar_values[0][6])
    print numpy.mean(tar_values[1][3]), numpy.mean(tar_values[1][4]), numpy.mean(tar_values[1][5]), numpy.mean(tar_values[1][6])
        
    X_n, Y_n = numpy.mgrid[0:1:201j, 0:1:200j]
    positions = numpy.vstack([X_n.ravel(), Y_n.ravel()])
    v_1 = numpy.vstack([tar_values[0][3], tar_values[0][4]])
    v_2 = numpy.vstack([tar_values[1][3], tar_values[1][4]])
    k_1 = gaussian_kde(v_1)
    k_2 = gaussian_kde(v_2)
    Z_1 = numpy.reshape(k_1(positions).T, X_n.shape)
    Z_2 = numpy.reshape(k_2(positions).T, X_n.shape)
    Z_3 = Z_2 - Z_1
    
    rmsd_square = 0.0
    for i_7 in range(len(Z_3)):
        for i_8 in range(len(Z_3[i_7])):
            rmsd_square += (Z_3[i_7][i_8])*(Z_3[i_7][i_8])
    print numpy.sqrt(rmsd_square) / 40200.0

    plot.figure(figsize = [7, 7])
    #plot.contour(numpy.rot90(numpy.fliplr(Z_1)), cmap = "Blues", levels = [1.0, 3.0, 9.0], extent = [0, 1, 0, 1])
    plot.imshow(numpy.rot90(Z_3), cmap = "seismic", extent = [0, 1, 0, 1], vmin = min(numpy.amin(Z_3), -numpy.amax(Z_3)), vmax = max(-numpy.amin(Z_3), numpy.amax(Z_3)))
        
    plot.plot([0, 1], [0.25, 0.25], 'k--')
    plot.plot([0, 1], [0.35, 0.35], 'k--')
    plot.plot([0.25, 0.25], [0, 1], 'k--')
    plot.plot([0.35, 0.35], [0, 1], 'k--')
    plot.xlabel("Pre-phosphorylation")
    plot.ylabel("Post-phosphorylation")
    plot.colorbar(cax = plot.axes([0.9, 0.1, 0.05, 0.8]))
    plot.legend(loc = 4)
    plot.show()
    plot.clf()
"""

####

"""
for i_1 in numpy.arange(1,2,1):
    tar_values = [[[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []]]
    for i_2 in range(5):
        for i_3 in range(len(vlist[i_1][i_2])):
            numeric = []
            nterm = 0
            cterm = 29
            for i_4 in range(29):
                det_1 = 0
                for i_5 in range(24):
                    if vlist[i_1][i_2][i_3][i_4] == aa_list[i_5]:
                        numeric.append(i_5)
                        det_1 += 1
                if det_1 < 1:
                    numeric.append(20)
                    det_1 += 1
                if vlist[i_1][i_2][i_3][28] != "_" and vlist[i_1][i_2][i_3][i_4] == "_":
                    nterm += 1
                if vlist[i_1][i_2][i_3][0] != "_" and vlist[i_1][i_2][i_3][i_4] == "_":
                    cterm -= 1
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
            sub_v = [numpy.mean(neg_pre[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.mean(neg_post[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.mean(pos_list[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.vdot(hp_pre[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(hp_post[(14-calc_range):(15+calc_range)], hp_weight)/121.0]
            tar_values[i_2][0].append(sub_v[0])
            tar_values[i_2][1].append(sub_v[1])
            tar_values[i_2][2].append(sub_v[2])
            tar_values[i_2][3].append(sub_v[2]-sub_v[0])
            tar_values[i_2][4].append(sub_v[2]-sub_v[1])
            tar_values[i_2][5].append(sub_v[3])
            tar_values[i_2][6].append(sub_v[4])

    #print numpy.mean(tar_values[1][3]), numpy.mean(tar_values[1][4]), numpy.mean(tar_values[1][5]), numpy.mean(tar_values[1][6])
        
    X_n, Y_n = numpy.mgrid[0:1:201j, 0:1:201j]
    positions = numpy.vstack([X_n.ravel(), Y_n.ravel()])
    v_1 = numpy.vstack([tar_values[0][3], tar_values[0][4]])
    v_2 = numpy.vstack([tar_values[1][3], tar_values[1][4]])
    v_3 = numpy.vstack([tar_values[2][3], tar_values[2][4]])
    v_4 = numpy.vstack([tar_values[3][3], tar_values[3][4]])
    v_5 = numpy.vstack([tar_values[4][3], tar_values[4][4]])
    k_1 = gaussian_kde(v_1)
    k_2 = gaussian_kde(v_2)
    k_3 = gaussian_kde(v_3)
    k_4 = gaussian_kde(v_4)
    k_5 = gaussian_kde(v_5)
    Z_x = [numpy.reshape(k_1(positions).T, X_n.shape), numpy.reshape(k_2(positions).T, X_n.shape), numpy.reshape(k_3(positions).T, X_n.shape), numpy.reshape(k_4(positions).T, X_n.shape), numpy.reshape(k_5(positions).T, X_n.shape)]
    
    for i_7 in range(len(target_pair)):
        Z_t = Z_x[target_pair[i_7][1]] - Z_x[target_pair[i_7][0]]
        rmsd_square = 0.0
        
        for i_9 in range(len(Z_t)):
            for i_10 in range(len(Z_t[i_9])):
                rmsd_square += (Z_t[i_9][i_10])*(Z_t[i_9][i_10])
        #print i_1, target_pair[i_7], numpy.fabs(numpy.sqrt(rmsd_square)) / 40000.0
        print i_1, target_pair[i_7], rmsd_square / 40000.0
        
        plot.figure(figsize = [7, 7])
        plot.imshow(numpy.rot90(Z_t), cmap = "bwr", extent = [0, 1, 0, 1], vmin = -16, vmax = 16)

        plot.plot([0, 1], [0.25, 0.25], 'k--')
        plot.plot([0, 1], [0.35, 0.35], 'k--')
        plot.plot([0.25, 0.25], [0, 1], 'k--')
        plot.plot([0.35, 0.35], [0, 1], 'k--')
        plot.xlabel("Pre-phosphorylation")
        plot.ylabel("Post-phosphorylation")
        plot.colorbar(cax = plot.axes([0.9, 0.1, 0.05, 0.8]))
        plot.legend(loc = 4)
        plot.show()
"""

####

for i_1 in numpy.arange(0,1,1):
    tar_values = [[[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], []]]
    for i_2 in range(5):
        for i_3 in range(len(vlist[i_1][i_2])):
            numeric = []
            nterm = 0
            cterm = 29
            for i_4 in range(29):
                det_1 = 0
                for i_5 in range(24):
                    if vlist[i_1][i_2][i_3][i_4] == aa_list[i_5]:
                        numeric.append(i_5)
                        det_1 += 1
                if det_1 < 1:
                    numeric.append(20)
                    det_1 += 1
                if vlist[i_1][i_2][i_3][28] != "_" and vlist[i_1][i_2][i_3][i_4] == "_":
                    nterm += 1
                if vlist[i_1][i_2][i_3][0] != "_" and vlist[i_1][i_2][i_3][i_4] == "_":
                    cterm -= 1
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
            #sub_v = [numpy.vdot(neg_pre[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(neg_post[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(pos_list[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(hp_pre[(14-calc_range):(15+calc_range)], hp_weight)/121.0, numpy.vdot(hp_post[(14-calc_range):(15+calc_range)], hp_weight)/121.0]
            sub_v = [numpy.mean(neg_pre[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.mean(neg_post[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.mean(pos_list[max(nterm, 14-calc_range):min(cterm, 15+calc_range)]), numpy.vdot(hp_pre[(14-calc_range):(15+calc_range)], hp_weight)/((float(sys.argv[1])+1)*(float(sys.argv[1])+1)), numpy.vdot(hp_post[(14-calc_range):(15+calc_range)], hp_weight)/((float(sys.argv[1])+1)*(float(sys.argv[1])+1))]
            tar_values[i_2][0].append(-sub_v[0])
            tar_values[i_2][1].append(-sub_v[1])
            tar_values[i_2][2].append(sub_v[2])
            tar_values[i_2][3].append(sub_v[2]-sub_v[0])
            tar_values[i_2][4].append(sub_v[2]-sub_v[1])
            #tar_values[i_2][3].append(numpy.fabs(sub_v[2]+sub_v[0]))
            #tar_values[i_2][4].append(numpy.fabs(sub_v[2]+sub_v[1]))
            tar_values[i_2][5].append(sub_v[3])
            tar_values[i_2][6].append(sub_v[4])
            #tar_values[i_2][5].append(sub_v[2]+sub_v[0])
            #tar_values[i_2][6].append(sub_v[2]+sub_v[1])

    #print numpy.mean(tar_values[0][3]), numpy.mean(tar_values[0][4]), numpy.mean(tar_values[0][5]), numpy.mean(tar_values[0][6])
    plot_titles = ["Negative charge - before phosphorylation", "Negative charge - after phosphorylation", "Positive charge", "Total charge - before phosphorylation", "Total charge - after phosphorylation", "Charge difference - before phosphorylation", "Charge difference - after phosphorylation"]
    plot_xlabels = ["S-nP", "T-nP", "Y", "S-P", "T-P"]
    
    fig, ax = plot.subplots(figsize = [8, 5])
    prop = fm.FontProperties(fname='/Users/cmdb/Desktop/phospho_analysis/190422/calibri/Calibri.ttf')
    rcParams.update({'font.size': 15})
    
    #ax.set_ylabel('Total charged AA fraction', fontproperties=prop, fontsize = 15) #
    ax.set_ylabel('PII propensity', fontproperties=prop, fontsize = 15) #
    ax.set_xlim([-0.5, 11.5])
    ax.set_ylim([0, 1])
    
    print kruskal(tar_values[0][3], tar_values[0][4])
    print kruskal(tar_values[1][3], tar_values[1][4])
    print kruskal(tar_values[2][3], tar_values[2][4])
    print kruskal(tar_values[3][3], tar_values[3][4])
    print kruskal(tar_values[4][3], tar_values[4][4])
    
    print kruskal(tar_values[0][3], tar_values[3][3])
    print kruskal(tar_values[0][4], tar_values[3][4])
    print kruskal(tar_values[1][3], tar_values[4][3])
    print kruskal(tar_values[1][4], tar_values[4][4])
    
    for i_7 in numpy.arange(5,7,1):
        tar_array = numpy.asarray([tar_values[0][i_7], tar_values[1][i_7], tar_values[2][i_7], tar_values[3][i_7], tar_values[4][i_7]])
        quart_1 = [numpy.percentile(tar_values[0][i_7], 25), numpy.percentile(tar_values[1][i_7], 25), numpy.percentile(tar_values[2][i_7], 25), numpy.percentile(tar_values[3][i_7], 25), numpy.percentile(tar_values[4][i_7], 25)]
        quart_3 = [numpy.percentile(tar_values[0][i_7], 75), numpy.percentile(tar_values[1][i_7], 75), numpy.percentile(tar_values[2][i_7], 75), numpy.percentile(tar_values[3][i_7], 75), numpy.percentile(tar_values[4][i_7], 75)]
        print numpy.mean(tar_values[0][i_7]), numpy.mean(tar_values[1][i_7]), numpy.mean(tar_values[2][i_7]), numpy.mean(tar_values[3][i_7]), numpy.mean(tar_values[4][i_7])
        print numpy.median(tar_values[0][i_7]), numpy.median(tar_values[1][i_7]), numpy.median(tar_values[2][i_7]), numpy.median(tar_values[3][i_7]), numpy.median(tar_values[4][i_7])
        
        if i_7 == 5:
            parts = ax.violinplot(tar_array, positions = numpy.arange(0, 12.5, 2.5), points = 500, showmeans = True, widths = 0.5)
            for pc in parts['bodies']:
                pc.set_facecolor('#C0C0C0')
                pc.set_edgecolor('#808080')
                pc.set_alpha(1)
            parts['cmeans'].set_color('#404040')
            parts['cmeans'].set_linewidth(3.0)
            parts['cmins'].set_color('#808080')
            parts['cmaxes'].set_color('#808080')
            parts['cbars'].set_color('#808080')
            ax.vlines(numpy.arange(0, 12.5, 2.5), quart_1, quart_3, color = '#404040', linestyle = '-', lw=3)
        elif i_7 == 6:
            parts = ax.violinplot(tar_array, positions = numpy.arange(1, 13.5, 2.5), points = 500, showmeans = True, widths = 0.5)
            for pc in parts['bodies']:
                pc.set_facecolor('#C0C0F0') #
                pc.set_edgecolor('#4040C0') #
                pc.set_alpha(1)
            parts['cmeans'].set_color('#0000B0') #
            parts['cmeans'].set_linewidth(3.0)
            parts['cmins'].set_color('#4040C0')
            parts['cmaxes'].set_color('#4040C0')
            parts['cbars'].set_color('#4040C0')
            ax.vlines(numpy.arange(1, 13.5, 2.5), quart_1, quart_3, color = '#0000B0', linestyle = '-', lw=3)
            
    #ax.plot([-1, 12.5], [0.25, 0.25], 'k--', linewidth = 1)
    #ax.plot([-1, 12.5], [0.35, 0.35], 'k--', linewidth = 1)
    ax.set_xticks(numpy.arange(0.5, 13, 2.5))
    ax.set_xticklabels(plot_xlabels)
    plot.xticks(fontproperties=prop, fontsize = 15)
    plot.yticks(fontproperties=prop, fontsize = 15)
    
    plot.show()

        