from __future__ import division

import numpy as np
from sklearn.ensemble import RandomForestClassifier
import inspect
import sys
import urllib2
import os
import re


subfldr = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"TCSS")))
sys.path.append(subfldr)
script_dir = subfldr.rsplit(os.sep, 1)[0] + os.sep
print script_dir

from main import load_semantic_similarity, calculate_semantic_similarity


class GOSim(object):
	objs = ""

	def __init__(self, taxid, gene_anno_F= script_dir + "go-basic.obo"):
		self.name = ["Sim_CC", "Sim_BP", "Sim_MF"]
		if GOSim.objs == "":
			GOSim.objs = load_semantic_similarity( gene_anno_F, taxid, "C:2.4,P:3.5,F:3.3", "")

	def calculateScore(self, edge):
		a, b = edge.split("\t")
		out = []
		domain_def = {'C': 'Cellular Component', 'P': 'Biological Process', 'F': 'Molecular Function'}
		for domain in domain_def:
			score = GOSim.objs[domain]._semantic_similarity(a, b)[0]
			if score is None: score = 0
			out.append(score)
		return out


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
		elm_FH.close()

	def countELMs(self, prot):
		counts = 0
		for elm in self.elms:
			hits = len(elm.findall(self.get_seq_prot(prot)))
			counts += hits
		self.prot2counts[prot] = counts

	def calculateScore(self, edge):
		protA, protB = edge.split("\t")
		for p in [protA, protB]:
			if p not in self.prot2counts: self.countELMs(p)
		return self.prot2counts[protA], self.prot2counts[protB]

class STRING:
	def __init__(self, taxid):
		self.name = ["Score", "Fusion", "Phylogenetic_profiling", "Homology", "Co-expression", "Experiment", "Database", "Textmining"]
		self.score_map = {"score" : 0 , "fscore" : 1, "pscore" : 2, "hscore" : 3, "ascore" : 4, "escore" : 5, "dscore" : 6, "tscore" : 7}
		self.taxid = taxid
		self.get_score = "http://string-db.org/api/psi-mi-tab/interactionsList?identifiers=%s%%0D%s&species=%s&limit=0"
		self.get_degree = "http://string-db.org/api/tsv/interactorsList?identifiers=%s%%0D%s&species=%s"

	def calculateScore(self, edge):
		a, b = edge.split("\t")
		out = [0] * len(self.score_map.keys())
 		scoreFH = urllib2.urlopen(self.get_score % (a,b, self.taxid))
		for line in scoreFH:
			line = line.rstrip()
			if line == "": continue
			scores = line.rsplit("\t",1)[1]
			for score in scores.split("|"):
				name, score = score.split(":")
				out[self.score_map[name]] = float(score)
		scoreFH.close()

		degree = 0
		degreeFH = urllib2.urlopen(self.get_degree % (a,b, self.taxid))
		degreeFH.readline()
		for _ in degreeFH:
			degree += 1
		degreeFH.close()
		out.append(degree-2)
		return out





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

class Degree:

		def __init__(self, networkfile = script_dir + "intact_biogrid.txt"):
			self.node_degree = {}
			netFH = open(networkfile)
			for line in netFH:
				line = line.rstrip()
				for node_id in line.split("\t"):
					if node_id not in self.node_degree:self.node_degree[node_id] = 0
					self.node_degree[node_id] += 1
			netFH.close()

		def calculateScore(self, edge):
			out = -2
			for node in edge.split("\t"):
				if node in self.node_degree: out += self.node_degree[node]
			return max(out, 0)



class CLF_Wrapper:
	def __init__(self):
		self.clf = RandomForestClassifier(n_estimators=100)
		self.fit()

	# Get reference data from online repository
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
	test = Degree()
	print test.calculateScore("O95622\tO60503")
	mode, infile = sys.argv[1:]
	if mode == "-test":
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