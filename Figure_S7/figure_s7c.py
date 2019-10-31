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

s_list = []
#s_list.extend(tar_list[0])
#s_list.extend(tar_list[1])
#s_list.extend(tar_list[2])
#s_list.extend(tar_list[3])
s_list.extend(tar_list[4])
tar_array = numpy.asarray(s_list)

plot.figure(figsize = [6, 6])
prop = fm.FontProperties(fname=file_path+'/External_Scripts/calibri/Calibri.ttf')
rcParams.update({'font.size': 15})

fpr_1, tpr_1, thres_1 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 1])
fpr_2, tpr_2, thres_2 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 2])
fpr_3, tpr_3, thres_3 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 3])
fpr_4, tpr_4, thres_4 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 4])
#fpr_5, tpr_5, thres_5 = metrics.roc_curve(tar_array[:, 0], tar_array[:, 5])

plot.plot([-0.1, 1.1], [-0.1, 1.1], 'k--')

plot.plot(fpr_4, tpr_4, linestyle = '-', color = '#000000', label = "PHOSforUS // AUC = %0.3f" % metrics.auc(fpr_4, tpr_4))
plot.plot(fpr_1, tpr_1, linestyle = '--', color = '#4040C0', label = "Horizontal features only // AUC = %0.3f" % metrics.auc(fpr_1, tpr_1))
plot.plot(fpr_2, tpr_2, linestyle = '--', color = '#C04040', label = "Vertical features only // AUC = %0.3f" % metrics.auc(fpr_2, tpr_2))
#plot.plot(fpr_3, tpr_3, linestyle = '--', color = '#808080', label = "PSWM // AUC = %0.3f" % metrics.auc(fpr_3, tpr_3))

plot.xlim([0, 1])
plot.ylim([0, 1])

plot.xlabel("False positive rate", fontproperties=prop, fontsize = 15)
plot.ylabel("True positive rate", fontproperties=prop, fontsize = 15)
plot.xticks(fontproperties=prop, fontsize = 15)
plot.yticks(fontproperties=prop, fontsize = 15)
plot.legend(loc=4, prop = prop, fontsize = '15')

plot.savefig("reproduced_figure_s7c.png")
