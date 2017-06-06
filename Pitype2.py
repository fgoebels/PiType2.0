from __future__ import division

import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.model_selection import  cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFECV

from sklearn import svm
from sklearn import metrics

# The following imported libraries are for intergrating GeneMANIA data as Functional Evidence
import sys
import urllib2
import os
import re

class ELM:
	def __init__(self):
		self.name = ["ELMA", "ELMB"]
		self.prot2seq = {}
		self.prot2counts = {}
		self.elms = []
		self.loadELM()

	def get_seq_prot(self, prot):
		if prot not in self.prot2seq: self.load_seq_prot(prot)
		return self.prot2seq[prot]

	def load_seq_prot(self, prot):
		fasta_url = "http://www.uniprot.org/uniprot/%s.fasta" % prot
		fasta_url_FH = urllib2.urlopen(fasta_url)
		fasta_url_FH.readline()
		seq = ""
		for line in fasta_url_FH.readline():
				line = line.rstrip()
				seq += line
		self.prot2seq[prot] = seq
		fasta_url_FH.close()

	def loadELM(self):
		elm_FH = urllib2.urlopen("http://elm.eu.org/elms/elms_index.tsv")
		for line in elm_FH:
			if line.startswith("#") or line.startswith("\"Accession\""): continue
			_, _, _, _, regex, _, _, _ = line.split("\t")
			if regex == "": continue
			regex = regex = re.sub("\"", "", regex)
			self.elms.append(re.compile(regex))

	def countELMs(self, prot):
		counts = 0
		for elm in self.elms:
			hits = len(elm.findall(self.get_seq_prot(prot)))
			counts += hits
		self.prot2counts[prot] = counts

	def getscore(self, edge):
		protA, protB = edge.split("\t")
		for p in [protA, protB]:
			if p not in self.prot2counts: self.countELMs(p)
		return self.prot2counts[protA], self.prot2counts[protB]

class Pitype:
	def __init__(self):
		return None


class CLF_Wrapper:
	def __init__(self, data, targets, forest=False):
		self.clf = RandomForestClassifier(n_estimators=100)

	def readData(self, file = ):


	def kFoldCV(self, folds=10):
		folds = StratifiedKFold(self.targets, folds)
		return cross_val_predict(self.clf, self.data, self.targets, cv=folds)

	def getValScores(self, folds=10):
		#		return cross_validation.cross_val_score(self.clf, self.data, self.targets, cv=10, scoring='f1')
		preds = self.kFoldCV(folds)
		precision = metrics.precision_score(self.targets, preds, average=None)[1]
		recall = metrics.recall_score(self.targets, preds, average=None)[1]
		fmeasure = metrics.f1_score(self.targets, preds, average=None)[1]
		auc_pr = average_precision_score(self.targets, preds)
		auc_roc = roc_auc_score(self.targets, preds)
		return [precision, recall, fmeasure, auc_pr, auc_roc]

	def getPRcurve(self, folds=10):
		all_probas = []
		all_targets = []
		for train, test in StratifiedKFold(self.targets, folds):
			probas = self.clf.fit(self.data[train], self.targets[train]).predict_proba(self.data[test])
			all_probas.extend(probas[:, 1])  # make sure that 1 is positive class in binarizied class vector
			all_targets.extend(self.targets[test])
		return precision_recall_curve(all_targets, all_probas)

	#		return roc_curve(all_targets, all_probas)

	def predict(self, toPred):
		return self.clf.predict_proba(toPred)

def main():
	thisELM = ELM()
	for edge in ["P04637\tP02340", "Q2XVY7\tQ29537"]:
		print thisELM.getscore(edge)
	return None

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass