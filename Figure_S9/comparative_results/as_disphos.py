#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.5
sv = []

for i_1 in range(5):
	fname = "disphos/dp_"+sys.argv[1]+"_"+str(i_1)+".txt"
	vlist = []
	plist = []
	slist = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			if int(ssl[0]) % 58 == 15:
				#vlist.append([int(ssl[0]), float(ssl[2])])
				vlist.append(float(ssl[2]))
				if float(ssl[2]) >= p_thres:
					plist.append(1)
				else:
					plist.append(0)
				slist.append(1)
			elif int(ssl[0]) % 58 == 44:
				#vlist.append([int(ssl[0]), float(ssl[2])])
				vlist.append(float(ssl[2]))
				if float(ssl[2]) >= p_thres:
					plist.append(1)
				else:
					plist.append(0)
				slist.append(0)
	#print len(vlist), len(plist), len(slist)
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

Sensitivity: 0.8333681618927521 0.04263939407470373
Specificity: 0.6615194384638274 0.018695137708845642
Precision  : 0.5439999999999999 0.035552777669262355
Neg-Pred-V : 0.89 0.03224903099319419
Accuracy   : 0.717 0.023151673805580433
F1-score   : 0.6574289016753675 0.03085600693412913
Matthews-CC: 0.4633727252262423 0.050052049917569265
AUROC      : 0.82301 0.028007166939910212

Sensitivity: 0.6326851927219843 0.03281438819963712
Specificity: 0.5554854165308327 0.018141128721031685
Precision  : 0.37 0.04147288270665544
Neg-Pred-V : 0.7859999999999999 0.020591260281973982
Accuracy   : 0.5780000000000001 0.02315167380558045
F1-score   : 0.4662286080576516 0.03934406743237973
Matthews-CC: 0.17127850011793802 0.04797704336199071
AUROC      : 0.6281100000000001 0.018086055401883544

Sensitivity: 0.6534310134310134 0.04956065464148873
Specificity: 0.5692032926336995 0.016101608967355963
Precision  : 0.41200000000000003 0.023151673805580447
Neg-Pred-V : 0.778 0.04791659420284376
Accuracy   : 0.595 0.024289915602982222
F1-score   : 0.504334451870477 0.022946554998832244
Matthews-CC: 0.2056000597757302 0.05595529686592245
AUROC      : 0.65703 0.024030139408667602

Sensitivity: 0.6919978886186491 0.06123957198919617
Specificity: 0.6891057893888084 0.044422361375764684
Precision  : 0.692 0.03599999999999999
Neg-Pred-V : 0.688 0.07573638491504597
Accuracy   : 0.6900000000000001 0.05244044240850758
F1-score   : 0.6914888842134792 0.04655401978719055
Matthews-CC: 0.3805509253568441 0.10493444872503459
AUROC      : 0.75849 0.0456587932385428

Sensitivity: 0.5683559619346622 0.01734427117620316
Specificity: 0.6441425106584062 0.05785239357320974
Precision  : 0.7619999999999999 0.054184868736576276
Neg-Pred-V : 0.42200000000000004 0.031240998703626618
Accuracy   : 0.592 0.026191601707417588
F1-score   : 0.6506316621781496 0.028274138475278533
Matthews-CC: 0.19760480021524862 0.06183326686251585
AUROC      : 0.6630299999999999 0.02285630328815227

"""

