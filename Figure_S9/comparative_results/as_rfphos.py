#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.5
sv = []

for i_1 in range(5):
	fname = "rfphos/rf_"+sys.argv[1]+"_"+str(i_1)+".txt"
	vlist = []
	plist = []
	slist = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			if len(ssl) > 3:
				if int(ssl[2]) % 58 == 15:
					vlist.append(float(ssl[3]))
					if float(ssl[3]) > p_thres:
						plist.append(1)
					else:
						plist.append(0)
					slist.append(1)
				elif int(ssl[2]) % 58 == 44:
					vlist.append(float(ssl[3]))
					if float(ssl[3]) > p_thres:
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

Sensitivity: 0.7910948319415139 0.017875952165621784
Specificity: 0.5900654198122552 0.01630593285356518
Precision  : 0.372 0.04445222154178575
Neg-Pred-V : 0.9020000000000001 0.011661903789690613
Accuracy   : 0.6369999999999999 0.019390719429665332
F1-score   : 0.5047148702840509 0.04096591628669832
Matthews-CC: 0.32294449771616 0.034386963671573324
AUROC      : 0.74387 0.025895030411258447

Sensitivity: 0.7112228338022353 0.045763907698040715
Specificity: 0.5750452491516848 0.01888602324016953
Precision  : 0.372 0.049558046773455475
Neg-Pred-V : 0.8479999999999999 0.03187475490101845
Accuracy   : 0.61 0.0244948974278318
F1-score   : 0.48663067644131236 0.04519549668253539
Matthews-CC: 0.2507196916360882 0.05269444221825179
AUROC      : 0.67351 0.03569695225085755

Sensitivity: 0.6233296456629789 0.0463969397741376
Specificity: 0.5762234539989893 0.028701460331572923
Precision  : 0.476 0.04223742416388577
Neg-Pred-V : 0.7120000000000001 0.03969886648255841
Accuracy   : 0.5940000000000001 0.034985711369071804
F1-score   : 0.5393518944481221 0.04177752411391939
Matthews-CC: 0.19367509108561287 0.07206388878133566
AUROC      : 0.63655 0.06700962617415497

Sensitivity: 0.5744101875803761 0.017476509927780125
Specificity: 0.6981501635262485 0.040449081291010325
Precision  : 0.836 0.024166091947189133
Neg-Pred-V : 0.38 0.03687817782917154
Accuracy   : 0.608 0.024207436873820428
F1-score   : 0.6808264376929741 0.01827236030985913
Matthews-CC: 0.24258380564725862 0.052375813302545864
AUROC      : 0.67445 0.03227266645320779

Sensitivity: 0.5528469065943181 0.02153157487186429
Specificity: 0.7011557231588287 0.08040488458875983
Precision  : 0.8699999999999999 0.04774934554525329
Neg-Pred-V : 0.296 0.048414873747640814
Accuracy   : 0.583 0.03280243893371344
F1-score   : 0.6757643258123229 0.02716836201197986
Matthews-CC: 0.20482114237040885 0.07884475231757049
AUROC      : 0.6442 0.03742290742312793

"""

