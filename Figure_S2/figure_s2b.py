#!/usr/bin/env python

import matplotlib.pyplot as plot
import numpy

vlabel = []
vlist = []
pcount = 0

with open("pre_calculated_data/log_190403.txt") as tar_file:
    for line in tar_file:
        ssl = line.rstrip('\n').split()
        if pcount < 220:
            if pcount < 120:
                vlabel.append(ssl[0]+ssl[1])
                vlist.append([float(ssl[2]), -float(ssl[3]), float(ssl[4]), float(ssl[5])])
            elif pcount >= 120:
                vlabel.append("+"+ssl[0]+ssl[1])
                vlist.append([float(ssl[2]), -float(ssl[3]), float(ssl[4]), float(ssl[5])])
        else:
            vlist[pcount-220].extend([float(ssl[2]), -float(ssl[3]), float(ssl[4]), float(ssl[5])])
        pcount += 1

varray = numpy.asarray(vlist)

v_1 = varray[varray[:,1].argsort()[::-1]]
v_2 = numpy.asarray(v_1)

xlabel = []
for i_1 in range(15):
    for i_2 in range(220):
        if vlist[i_2][1] == v_2[i_1][1]:
            xlabel.append(vlabel[i_2])
            break

x_index = numpy.arange(0.0, 15.0, 1.0)

plot.figure(figsize = [12, 7])
plot.bar(x_index-0.2, v_2[:15, 1], 0.3, color = 'blue', yerr = v_2[:15, 3], label = "Phosphorylation sites")
plot.bar(x_index+0.2, v_2[:15, 5], 0.3, color = 'orange', yerr = v_2[:15, 7], label = "Non-phosphorylated sequences")

plot.legend(loc=1)
plot.xlim([-0.5, 14.5])
plot.xticks(x_index, xlabel)
plot.xlabel("Amino acid criteria")
plot.ylabel("-log(p-value)")
plot.savefig("reproduced_figure_s2b.png")
