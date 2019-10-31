#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.5
sv = []

for i_1 in range(5):
	fname = "phprrf/pr_"+sys.argv[1]+"_"+str(i_1)+".txt"
	raw_list = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			raw_list.append([int(ssl[1]), float(ssl[3])])

	vlist = []
	plist = []
	slist = []
	clist = []

	for i_2 in range(100):
		tar_p = i_2*58+15
		tar_n = i_2*58+44
		det_p = 0
		det_n = 0
		for i_3 in range(len(raw_list)):
			if raw_list[i_3][0] == tar_p:
				det_p += 1
				vlist.append(raw_list[i_3][1])
				plist.append(1)
				slist.append(1)
				clist.append(tar_p)
				break
		if det_p < 1:
			vlist.append(0.0)
			plist.append(0)
			slist.append(1)
			clist.append(tar_p)
		for i_4 in range(len(raw_list)):
			if raw_list[i_4][0] == tar_n:
				det_n += 1
				vlist.append(raw_list[i_3][1])
				plist.append(1)
				slist.append(0)
				clist.append(tar_n)
		if det_n < 1:
			vlist.append(0.0)
			plist.append(0)
			slist.append(0)
			clist.append(tar_n)

	print len(vlist), len(plist), len(slist), len(clist)
	print clist[0:10]
	confmat = metrics.confusion_matrix(plist, slist).ravel()
	tpr = float(confmat[3])/float(confmat[2]+confmat[3])
	tnr = float(confmat[0])/float(confmat[0]+confmat[1])
	ppv = float(confmat[3])/float(confmat[1]+confmat[3])
	npv = float(confmat[0])/float(confmat[0]+confmat[2])
	acc = metrics.accuracy_score(slist, plist)
	f1 = metrics.f1_score(slist, plist)
	mcc = metrics.matthews_corrcoef(slist, plist)
	auroc = metrics.roc_auc_score(slist, vlist)
	sv.append([tpr, tnr, ppv, npv, acc, f1, mcc, auroc])

st = numpy.transpose(sv)

print "Sensitivity:", numpy.mean(st[0]), numpy.std(st[0])
print "Specificity:", numpy.mean(st[1]), numpy.std(st[1])
print "Precision  :", numpy.mean(st[2]), numpy.std(st[2])
print "Neg-Pred-V :", numpy.mean(st[3]), numpy.std(st[3])
print "Accuracy   :", numpy.mean(st[4]), numpy.std(st[4])
print "F1-score   :", numpy.mean(st[5]), numpy.std(st[5])
print "Matthews-CC:", numpy.mean(st[6]), numpy.std(st[6])
print "AUROC      :", numpy.mean(st[7]), numpy.std(st[7])

"""

Sensitivity: 0.8577540360873694 0.05270954140119425
Specificity: 0.72109452181016 0.03799846027098853
Precision  : 0.654 0.05953150426454885
Neg-Pred-V : 0.89 0.046043457732885346
Accuracy   : 0.772 0.041060930335295655
F1-score   : 0.7405944870372901 0.05000984369104305
Matthews-CC: 0.5610349204219425 0.08138953096134358
AUROC      : 0.77431 0.03420009356712349

Sensitivity: 0.7013851497108148 0.0254504449159069
Specificity: 0.6414951075164724 0.022271269631057675
Precision  : 0.5780000000000001 0.03370459909270543
Neg-Pred-V : 0.754 0.023323807579381222
Accuracy   : 0.666 0.023323807579381184
F1-score   : 0.6334648959629268 0.02842218472557252
Matthews-CC: 0.3373891926662006 0.0464645240774514
AUROC      : 0.66585 0.028836487303414815

Sensitivity: 0.5964242033831628 0.02330489250991057
Specificity: 0.6589304321409585 0.035212321878523474
Precision  : 0.744 0.021540659228538036
Neg-Pred-V : 0.496 0.03555277766926237
Accuracy   : 0.6199999999999999 0.028106938645110418
F1-score   : 0.6620588054590886 0.022749677275548757
Matthews-CC: 0.24755634782313624 0.05733677047835797
AUROC      : 0.6358499999999999 0.022187631689750017

Sensitivity: 0.6022021080798764 0.020690555950641905
Specificity: 0.8821291495321534 0.016851272987048405
Precision  : 0.95 0.012649110640673528
Neg-Pred-V : 0.37 0.06164414002968975
Accuracy   : 0.66 0.025690465157330252
F1-score   : 0.7367702736761025 0.012723212843422508
Matthews-CC: 0.3928395190888089 0.039483921864287716
AUROC      : 0.66763 0.02665199054479799

Sensitivity: 0.5535522026640909 0.010168901792944917
Specificity: 0.8445038621509209 0.03927128344338441
Precision  : 0.9559999999999998 0.01854723699099139
Neg-Pred-V : 0.22800000000000004 0.04308131845707604
Accuracy   : 0.592 0.015033296378372921
F1-score   : 0.7009412952610081 0.006467076739345044
Matthews-CC: 0.26936559739733607 0.02576223323899824
AUROC      : 0.5872399999999999 0.022632971523863142

"""

