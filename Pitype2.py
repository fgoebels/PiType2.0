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

	def __init__(self, fs):

		self.clf = CLF_Wrapper()
		self.fs = self.get_fs_comb(fs)
		return None

	def get_fs_comb(self, comb_string):
		# Create feature combination
		scores = [ELM()]
		this_scores = []
		for i, feature_selection in enumerate(comb_string):
			if feature_selection == "1": this_scores.append(scores[i])
		return this_scores



class CLF_Wrapper:
	def __init__(self):
		self.clf = RandomForestClassifier(n_estimators=100)
		self.fit()

	def fit(self, data_loc = "https://raw.githubusercontent.com/fgoebels/PiType2.0/master/ob_nob.dat.txt"):
		data_FH = urllib2.urlopen(data_loc)
		data_FH.readline()
		targets = []
		data = []
		for line in data_FH:
			line = line.rstrip()
			line = line.split("\t")
			scores = line[2:-1]
			data.append(map(float, scores))
			label = line[-1]
			if label == "nob":
				targets.append(0)
			else:
				targets.append(1)
		data_FH.close()
		targets = np.array(targets)
		data = np.array(data)
		self.clf.fit(data, targets)

	def predict_proba(self, toPred):
		probas = self.clf.predict_proba(toPred)
		return probas[:,1]

	def predict(self, toPred):
		preds = self.clf.predict(toPred)
		return preds

def main():
	this_ELM = ELM()
	this_clf =  CLF_Wrapper()
	for edge in ["P04637\tP02340", "Q2XVY7\tQ29537"]:
		scores = np.array(this_ELM.getscore(edge)).reshape(1, -1)
		print this_clf.predict(scores)
		print this_clf.predict_proba(scores)
	return None

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass