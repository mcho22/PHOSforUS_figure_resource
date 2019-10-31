#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.5
sv = []

for i_1 in range(5):
	fname = "netphos/np_"+sys.argv[1]+"_"+str(i_1)+".txt"
	clist = [""]
	vlist = []
	plist = []
	slist = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			if len(ssl) > 1:
				if ssl[2] != clist[-1]:
					clist.append(ssl[2])
					if int(ssl[2]) % 58 == 15:
						vlist.append(float(ssl[5]))
						if float(ssl[5]) >= p_thres:
							plist.append(1)
						else:
							plist.append(0)
						slist.append(1)
					elif int(ssl[2]) % 58 == 44:
						vlist.append(float(ssl[5]))
						if float(ssl[5]) >= p_thres:
							plist.append(1)
						else:
							plist.append(0)
						slist.append(0)

	print len(vlist), len(plist), len(slist)
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

Sensitivity: 0.5779899823461163 0.018461548772474617
Specificity: 0.7285668276972623 0.04401547988196343
Precision  : 0.8619999999999999 0.027129319932501096
Neg-Pred-V : 0.37 0.03898717737923585
Accuracy   : 0.6159999999999999 0.025573423705088867
F1-score   : 0.6918467144876003 0.0195981341528459
Matthews-CC: 0.2665709433509945 0.055292351102725
AUROC      : 0.7171100000000001 0.043992867603737794

Sensitivity: 0.5281683784624962 0.01925348232676364
Specificity: 0.5345646074269262 0.022048022650328805
Precision  : 0.5980000000000001 0.00979795897113272
Neg-Pred-V : 0.46399999999999997 0.04586937976471887
Accuracy   : 0.531 0.02059126028197402
F1-score   : 0.5606287560978696 0.010012628989527584
Matthews-CC: 0.06236466765157643 0.041233813699962894
AUROC      : 0.53124 0.04163952929609076

Sensitivity: 0.6169332333083271 0.021491299940192336
Specificity: 0.5933270372192162 0.00970379201451961
Precision  : 0.55 0.03847076812334268
Neg-Pred-V : 0.6559999999999999 0.04841487374764082
Accuracy   : 0.603 0.011224972160321832
F1-score   : 0.5801753858153658 0.017003363692770664
Matthews-CC: 0.20810317712970786 0.02406022568872331
AUROC      : 0.6224700000000001 0.02159112317597212

Sensitivity: 0.5171750850292925 0.013129899589980818
Specificity: 0.751391223155929 0.18059526969972722
Precision  : 0.9720000000000001 0.017204650534085267
Neg-Pred-V : 0.092 0.03655133376499413
Accuracy   : 0.532 0.0242074368738204
F1-score   : 0.675090297158107 0.014589367758369078
Matthews-CC: 0.13008929386776727 0.0936401581327112
AUROC      : 0.63346 0.05346838692161941

Sensitivity: 0.5107478710596064 0.006781906064272952
Specificity: 0.646410899042478 0.09249814696726044
Precision  : 0.95 0.020976176963403006
Neg-Pred-V : 0.09 0.02
Accuracy   : 0.5199999999999999 0.012649110640673528
F1-score   : 0.6642868346643702 0.009871361110923182
Matthews-CC: 0.07912789910896707 0.05011307565647626
AUROC      : 0.59619 0.03318760913353051

"""

