#!/usr/bin/env python

#### Generic module import
import datetime
import matplotlib.pyplot as plot
import numpy
import random
import sys

value_list = [[0.882004, 0.837234, 0.850022, 0.832020, 0.805172], [0.858714, 0.906428, 0.829778, 0.883542, 0.812794], [0.822862, 0.791984, 0.860508, 0.793345, 0.746766], [0.779316, 0.809318, 0.768972, 0.834062, 0.737062], [0.738412, 0.734118, 0.727142, 0.739020, 0.821618]]

p_class = ["SP", "S-", "TP", "T-", "Y-"]

fig, axis = plot.subplots(figsize = (7, 7))
heatmap = axis.pcolor(value_list, cmap=plot.cm.OrRd)

plot.colorbar(heatmap)
plot.xlabel("Applied PSWM sets")
plot.ylabel("Target sequence subclass")
plot.xticks(numpy.arange(0.5, 5.5, 1.0), p_class)
plot.yticks(numpy.arange(0.5, 5.5, 1.0), p_class)

for i in range(5):
    for j in range(5):
        if value_list[i][j] >= 0.8:
            text = axis.text(float(j)+0.5, float(i)+0.5, str(value_list[i][j])[0:6], ha="center", va="center", color="w")
        else:
            text = axis.text(float(j)+0.5, float(i)+0.5, str(value_list[i][j])[0:6], ha="center", va="center", color="k")
plot.savefig("reproduced_figure_s2c.png")


