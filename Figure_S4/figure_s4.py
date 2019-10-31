#!/usr/bin/env python

import matplotlib.pyplot as plot
import matplotlib.font_manager as fm
from matplotlib import rcParams
import numpy
import sys

from scipy.optimize import curve_fit, differential_evolution
from scipy.stats import gaussian_kde, ttest_ind, f_oneway, kruskal

########

file_path = "/Users/cmdb/Desktop/phospho_analysis/191007_figurebase" # Change as appropriate

########

base_l = 2.16*numpy.sqrt(6)
base_aa = 29

def sigmoid(x, x0, k, s_a, s_b):
    y = s_a/(1.0+numpy.exp(-x0*(x-k)))+s_b
    return y

base_min = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.0))*base_l
base_max = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.95))*base_l
base_diff = base_max - base_min

b_20 = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.2))*base_l
b_40 = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.4))*base_l
b_60 = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.6))*base_l
b_80 = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-0.8))*base_l

#zarray = (numpy.asarray([[0.0, 0.0], [0.2, 0.0],[0.25, 0.01], [0.2727, 0.01], [0.3333, 0.023], [0.3333, 0.051], [0.428571, 0.140], [0.5, 0.160], [0.6, 0.287], [1.0, 0.609]]))
zarray = (numpy.asarray([[0.0, 0.0], [0.2, 0.0],[0.25, 0.01], [0.3333, 0.023], [0.5, 0.160], [0.6, 0.287], [1.0, 0.609]]))
popt, pcov = curve_fit(sigmoid, zarray[:, 0], zarray[:, 1])

plot_xrange = numpy.arange(20, 60.5, 0.5)
plot_yrange = numpy.arange(0.0, 40.5, 0.5)
plot_zval = []

for i_1 in plot_yrange:
    #base_charge = base_min + base_diff * max(0.0, sigmoid(i_1*0.01, *popt))
    base_charge = base_min + base_diff * sigmoid(i_1*0.01, *popt)
    plot_zsub = []
    for i_2 in plot_xrange:
        #b_val = ((1.24*i_2+0.904)*(0.00759*i_1*(float(base_l))+0.963)*2.49)*numpy.power(base_aa, 0.509)*numpy.sqrt(6)
        #plot_zsub.append(b_val)
        b_val = numpy.power(base_aa, 0.503-0.11*numpy.log(1.0-i_2*0.01))*base_l
        plot_zsub.append(base_charge+(base_max-base_charge)*((b_val-base_min)/base_diff))
    plot_zval.append(plot_zsub)

plot_zarray = numpy.asarray(plot_zval)

fig, ax = plot.subplots(figsize = [9, 8])
prop = fm.FontProperties(fname=file_path+'/External_Scripts/calibri/Calibri.ttf')
rcParams.update({'font.size': 15})

cs = ax.contourf(plot_xrange, plot_yrange, plot_zarray, numpy.arange(30, 40.1, 0.1), cmap=plot.get_cmap("RdBu_r"))
ax.contour(plot_xrange, plot_yrange, plot_zarray, numpy.arange(30, 40.5, 0.5), colors = 'w', linestyles='dotted')
cbar = fig.colorbar(cs)

x_add = 0.0
plot.plot([35.3-x_add, 36.0-x_add], [12.68, 16.21], 'ro', markersize = 5.0, alpha = 0.5)
plot.plot([43.0-x_add, 44.6-x_add], [8.89, 10.78], 'bo', markersize = 5.0, alpha = 0.5)
plot.plot([36.2-x_add, 35.4-x_add], [10.63, 16.07], 'ko', markersize = 5.0, alpha = 0.5)

plot.arrow(35.3-x_add, 12.68, 0.7, 3.53, color = 'r', width = 0.3, head_width = 1.2, head_length = 1.2, length_includes_head = True)
plot.arrow(43.0-x_add, 8.89, 1.6, 1.90, color = 'b', width = 0.3, head_width = 1.2, head_length = 1.2, length_includes_head = True)
plot.arrow(36.2-x_add, 10.63, -0.8, 5.44, color = 'k', width = 0.3, head_width = 1.2, head_length = 1.2, length_includes_head = True)

ax.set_xlim([25, 55])
ax.set_ylim([0, 30])
ax.set_xlabel("PII propensity (%)", fontproperties=prop, fontsize = 16)
ax.set_ylabel("Absolute net charge (%)", fontproperties=prop, fontsize = 16)
ax.tick_params(labelsize = 15)

plot.savefig("reproduced_figure_s4.png")

