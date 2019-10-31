#!/usr/bin/env python

from matplotlib import rcParams
import matplotlib.font_manager as fm
import matplotlib.pyplot as plot
import numpy

from sklearn import metrics

########

file_path = "/Users/cmdb/Desktop/phospho_analysis/191007_figurebase" # Change as appropriate

########

tar_list = [[], [], [], [], []]

tar_file = open("roc_recalc.txt")

for line in tar_file:
    ssl = line.rstrip('\n').split('\t')
    ssl_2 = [int(ssl[1]), float(ssl[2]), float(ssl[3]), float(ssl[4]), float(ssl[5])]
    for i_0 in range(4):
        if float(ssl[1+i_0]) > numpy.log(0.5):
            ssl_2.append(1)
        else:
            ssl_2.append(0)
    tar_list[int(ssl[0])].append(ssl_2)

tar_file.close()

plen = [len(tar_list[0])/10, len(tar_list[1])/10, len(tar_list[2])/10, len(tar_list[3])/10, len(tar_list[4])/10]
vlist = []
"""
for i_1 in range(10):
    sublist = []
    for i_2 in range(5): #
        for i_3 in range(plen[i_2]):
            sublist.append(tar_list[i_2][(i_1*plen[i_2])+i_3])
    subarray = numpy.asarray(sublist)
    #vsub = []
    #for i_4 in range(4):
        #tn, fp, fn, tp = metrics.confusion_matrix(subarray[:, 0], subarray[:, 5+i_4]).ravel()
        #vsub.append(float(tn)/float(tn+fp))
    #vlist.append(vsub)
    vsub = [metrics.roc_auc_score(subarray[:, 0], subarray[:, 1]), metrics.roc_auc_score(subarray[:, 0], subarray[:, 2]), metrics.roc_auc_score(subarray[:, 0], subarray[:, 3]), metrics.roc_auc_score(subarray[:, 0], subarray[:, 4])]
    vlist.append(vsub)
    print len(sublist), vsub

varray = numpy.asarray(vlist)
print numpy.mean(varray[:, 0]), numpy.std(varray[:, 0])
print numpy.mean(varray[:, 1]), numpy.std(varray[:, 1])
print numpy.mean(varray[:, 2]), numpy.std(varray[:, 2])
print numpy.mean(varray[:, 3]), numpy.std(varray[:, 3])
"""
s_list = []
s_list.extend(tar_list[0])
s_list.extend(tar_list[1])
s_list.extend(tar_list[2])
s_list.extend(tar_list[3])
s_list.extend(tar_list[4])
tar_array = numpy.asarray(s_list)

plot.figure(figsize = [6, 6])
prop = fm.FontProperties(fname=file_path+'/External_Scripts/calibri/Calibri.ttf')
rcParams.update({'font.size': 15})

#fpr_1, tpr_1, thres_1 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 1])
#fpr_2, tpr_2, thres_2 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 2])
#fpr_3, tpr_3, thres_3 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 3])
#fpr_4, tpr_4, thres_4 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 4])
#fpr_5, tpr_5, thres_5 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 5])
fpr_1, tpr_1, thres_1 = metrics.precision_recall_curve(tar_array[:, 0], tar_array[:, 1])
fpr_2, tpr_2, thres_2 = metrics.precision_recall_curve(tar_array[:, 0], tar_array[:, 2])
fpr_3, tpr_3, thres_3 = metrics.precision_recall_curve(tar_array[:, 0], tar_array[:, 3])
fpr_4, tpr_4, thres_4 = metrics.precision_recall_curve(tar_array[:, 0], tar_array[:, 4])
#fpr_5, tpr_5, thres_5 = metrics.precision_recall_curve(tar_array[:, 0], tar_array[:, 5])

plot.plot(tpr_4, fpr_4, linestyle = '-', color = '#000000', label = "PHOSforUS / All features")
plot.plot(tpr_1, fpr_1, linestyle = '--', color = '#4040C0', label = "Horizontal features only")
plot.plot(tpr_2, fpr_2, linestyle = '--', color = '#C04040', label = "Vertical features only")
plot.plot(tpr_3, fpr_3, linestyle = '--', color = '#808080', label = "Position-specific weight matrix")

plot.xlim([0, 1])
plot.ylim([0, 1])

plot.xlabel("Recall", fontproperties=prop, fontsize = 15)
plot.ylabel("Precision", fontproperties=prop, fontsize = 15)
plot.xticks(fontproperties=prop, fontsize = 15)
plot.yticks(fontproperties=prop, fontsize = 15)
plot.legend(loc=4, prop = prop, fontsize = '15')

plot.savefig("reproduced_figure_s8b.png")
