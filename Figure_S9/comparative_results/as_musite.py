#!/usr/bin/env python

import numpy
import sys

from sklearn import metrics

####

p_thres = 0.0
sv = []

for i_1 in range(5):
	fname = "musite/ms_"+sys.argv[1]+"_"+str(i_1)+".txt"
	vlist = []
	plist = []
	slist = []
	with open(fname) as p_file:
		for line in p_file:
			ssl = line.rstrip('\n').split()
			if len(ssl) > 3:
				if int(ssl[0]) % 58 == 15:
					vlist.append(float(ssl[3]))
					if float(ssl[3]) >= p_thres:
						plist.append(1)
					else:
						plist.append(0)
					slist.append(1)
				elif int(ssl[0]) % 58 == 44:
					vlist.append(float(ssl[3]))
					if float(ssl[3]) >= p_thres:
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

Sensitivity: 0.8041456139125396 0.04641294590297294
Specificity: 0.6174275187202088 0.01743585643320318
Precision  : 0.44800000000000006 0.03709447398198281
Neg-Pred-V : 0.89 0.030983866769659363
Accuracy   : 0.669 0.02416609194718914
F1-score   : 0.5745165441057469 0.03658702637557355
Matthews-CC: 0.37737915996108795 0.053510119528111286
AUROC      : 0.7834600000000002 0.036931049267520145

Sensitivity: 0.6849026891706371 0.042658800042838826
Specificity: 0.5681609835007824 0.01955518858357551
Precision  : 0.36599999999999994 0.054990908339470075
Neg-Pred-V : 0.8320000000000001 0.030594117081556692
Accuracy   : 0.599 0.026153393661244032
F1-score   : 0.475075468802941 0.053479007000535095
Matthews-CC: 0.22363731199405942 0.05519971254806012
AUROC      : 0.67413 0.018760612996381522

Sensitivity: 0.6063589816642623 0.03756402059554109
Specificity: 0.5951433640531385 0.031089423411448266
Precision  : 0.578 0.036551333764994115
Neg-Pred-V : 0.6220000000000001 0.05979966555090423
Accuracy   : 0.5999999999999999 0.03361547262794321
F1-score   : 0.5909977665182501 0.030217805534977396
Matthews-CC: 0.20074825393683593 0.06736615741907058
AUROC      : 0.65107 0.03890687085850005

Sensitivity: 0.5910412152908853 0.03246420916933007
Specificity: 0.7486271111750707 0.061148383655357054
Precision  : 0.868 0.04489988864128729
Neg-Pred-V : 0.394 0.09308061022576077
Accuracy   : 0.631 0.040669398815325546
F1-score   : 0.702042159294411 0.02499329431478168
Matthews-CC: 0.29766203881575437 0.08196732924752627
AUROC      : 0.71465 0.029883724667450654

Sensitivity: 0.5653984319182137 0.014212976559822009
Specificity: 0.6879672623766488 0.04360634125745295
Precision  : 0.8379999999999999 0.025612496949731386
Neg-Pred-V : 0.35600000000000004 0.019595917942265426
Accuracy   : 0.597 0.02135415650406264
F1-score   : 0.6752046105612184 0.018172543454667317
Matthews-CC: 0.22169177674720636 0.049543366011138894
AUROC      : 0.6354 0.010523450004632489

"""

