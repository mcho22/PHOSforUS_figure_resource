#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.5
sv = []

for i_1 in range(5):
	fname = "phosphosvm/ps_"+sys.argv[1]+"_"+str(i_1)+".txt"
	vlist = []
	plist = []
	slist = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			if int(ssl[1]) % 58 == 15:
				#vlist.append([int(ssl[0]), float(ssl[2])])
				vlist.append(float(ssl[3]))
				if float(ssl[3]) >= p_thres:
					plist.append(1)
				else:
					plist.append(0)
				slist.append(1)
			elif int(ssl[1]) % 58 == 44:
				#vlist.append([int(ssl[0]), float(ssl[2])])
				vlist.append(float(ssl[3]))
				if float(ssl[3]) >= p_thres:
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

Sensitivity: 0.8737514940991847 0.03978118680405952
Specificity: 0.598808015416411 0.008732779375797388
Precision  : 0.366 0.023323807579381198
Neg-Pred-V : 0.946 0.020591260281973982
Accuracy   : 0.6559999999999999 0.012806248474865708
F1-score   : 0.5151236572730747 0.023041529438436865
Matthews-CC: 0.3837914442935248 0.03160118781247798
AUROC      : 0.8135600000000001 0.025226006421944804

Sensitivity: 0.7794231513743709 0.04286424170787129
Specificity: 0.5632855409661668 0.010005862240141406
Precision  : 0.288 0.027129319932501072
Neg-Pred-V : 0.9180000000000001 0.019390719429665297
Accuracy   : 0.603 0.015684387141358135
F1-score   : 0.41985562761532724 0.03169216429531511
Matthews-CC: 0.26553390176390534 0.038929616413768416
AUROC      : 0.7200199999999999 0.02440931379617217

Sensitivity: 0.6049834731756929 0.00868263210585845
Specificity: 0.6375780645928747 0.012444501019466358
Precision  : 0.686 0.017435595774162656
Neg-Pred-V : 0.5519999999999999 0.015999999999999976
Accuracy   : 0.619 0.009695359714832666
F1-score   : 0.6428543962014946 0.010520436321616356
Matthews-CC: 0.24026872062389856 0.01970408242792951
AUROC      : 0.67874 0.020602096980647355

Sensitivity: 0.5288969367677284 0.012161360897352709
Specificity: 0.867531574033122 0.06052219236427052
Precision  : 0.9800000000000001 0.015491933384829681
Neg-Pred-V : 0.12599999999999997 0.047159304490206375
Accuracy   : 0.5529999999999999 0.021354156504062596
F1-score   : 0.6868779426544149 0.010461907839832212
Matthews-CC: 0.20144880366180756 0.05413038197740759
AUROC      : 0.70333 0.05436889368011824

Sensitivity: 0.5519871791408943 0.004808195604085746
Specificity: 0.8478730158730159 0.06427372806148438
Precision  : 0.9560000000000002 0.02870540018881463
Neg-Pred-V : 0.22400000000000003 0.026532998322843202
Accuracy   : 0.5899999999999999 0.008366600265340763
F1-score   : 0.6997219251336898 0.009211928056561177
Matthews-CC: 0.26751991140219017 0.03073135065744997
AUROC      : 0.69151 0.0298511205819815

"""

