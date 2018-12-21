import csv
import copy
import time
import numpy as np
import pandas as pd
import itertools
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt

#不省略顯示
np.set_printoptions(threshold=np.inf)

def convertIBMdata():
	IBM = open( 'data.ntrans_0.1.tlen_40.nitems_0.1.txt','r')
	f = open('graph_7.txt','w')
	lines=IBM.readlines()

	for line in lines:
		print(line.split()[1]+','+line.split()[2])
		f.write(line.split()[1]+','+line.split()[2]+'\n')
		pass
	IBM.close()
	f.close()

def convertRuledata():
	Rule = open( 'association_rules.txt','r')
	f = open('graph_8.txt','w')
	lines=Rule.readlines()

	for line in lines:
		print(line.split()[0]+','+line.split()[1])
		f.write(line.split()[0]+','+line.split()[1]+'\n')
		pass
	Rule.close()
	f.close()

def edgelist(fileName):
	reader = csv.reader( open(fileName, 'rt') )
	return [[item for item in row ] for row in reader] #return edgelist

def Nodenum(edgeList):
	nodeList = []
	for edge in edgeList:
		for node in edge:
			if node not in nodeList:
				nodeList.append(node)
	return len(nodeList)

def nx_Graph(edgeList):
	# create null direted graph G
	G = nx.DiGraph ()
	for edge in edgeList :
		G.add_edge( int(edge[0]), int(edge[1]) )
	return G

def nx_HITS(G):
	hub, auth = nx.hits(G, max_iter = 1000, tol = 0.000005)
	return hub, auth

def nx_PageRank(G,d):
	PR = nx.pagerank(G, alpha = d)
	return PR

def DisplayGraph(G):
	pos = nx.spring_layout(G)
	nx.draw(G,node_size =300,node_color='pink',style='dashed',width=2.0,edge_color='gray')
	plt.show()

def adjMetrix(edgeList):
	m = Nodenum(edgeList)
	matrix = np.zeros(shape=(m,m))
	for edge in edgeList:

		u = int(edge[0]) - 1
		v = int(edge[1]) - 1
		matrix[u][v] = 1
	return matrix

def HITS(adjmatrix):
	m,m = adjmatrix.shape

	#hub 值都先設1
	hub = np.ones(shape=(m,1))
	hub_prev = np.zeros(shape=(m,1))
	#print(hub,hub_prev)

	#authority = A^T * hub, hub = A * authority
	authority = np.dot(adjmatrix.transpose(),hub)
	hub = np.dot(adjmatrix,authority)

	# check converge
	while not Converge_checking(hub, hub_prev):
		hub_prev = hub

		authority = np.dot(adjmatrix.transpose(),hub)
		hub = np.dot(adjmatrix,authority)	

		# normalize
		hub = hub/hub.sum()
		authority = authority/authority.sum()

	return hub, authority

def PageRank(adjmatrix,d):
	m,m = adjmatrix.shape


	outdegree = adjmatrix.sum(axis=1)
	PR = np.zeros(shape=(m,1))
	PR_prev = np.zeros(shape=(m,1))
	PR_Div_outdeg = np.zeros(shape=(m,1))
	PR_Div_outdeg_sum = np.zeros(shape=(m,1))

	PR = [ 1/m for value in PR]


	while not Converge_checking(PR, PR_prev):
		PR_prev = PR


		for i in range(0,m):
			if outdegree[i]!=0:
				PR_Div_outdeg[i] = PR[i] / outdegree[i]

		for i in range(0,m):
			PR_Div_outdeg_sum_value = 0
			for j in range(0,m):
				if adjmatrix[j][i] == 1:
					PR_Div_outdeg_sum_value +=  PR_Div_outdeg[j]
			PR_Div_outdeg_sum[i] = PR_Div_outdeg_sum_value

		PR = d / m + (1-d) * PR_Div_outdeg_sum
	#normolized
	PR = PR/PR.sum()
	return PR

def SimRank(G,n,C=0.8,t=100):
	
	S = np.identity(n)
	I = np.identity(n)
	G = G/G.sum(0)
	i = 1

	# S(a,b) ifa=b S設1  
	for a in range(t):
		S = C * np.dot(np.dot(G.T,S),G) + (1-C) * I
	for j in range(n):
		S[j][j] = 1

	return S
	
def Converge_checking(vector,vector_pre):
	converge = False
	tolerance = 0.000005
	changes = 0
	m = len(vector)

	for i in range(0,m):
		changes += abs(vector[i] - vector_pre[i])

	if changes < tolerance:
		converge = True
	return converge

def DisplayResult(hub,authority,PR,SR):


	Result_dict = {	
	                'Hub': hub.tolist(),
	                'Authority': authority.tolist(),
	                'PageRank' :PR.tolist(),
	                'SimRank' :SR.tolist()
	}

	Result_df = pd.DataFrame(Result_dict)
	#印出全部矩陣
	pd.set_option('display.max_columns',None)
	pd.set_option('display.max_rows', None)
	pd.set_option('display.max_colwidth',200)
	pd.set_option('precision', 4)
	#不省略顯示
	np.set_printoptions(threshold=np.inf)



	print(Result_df[['Hub','Authority']].sort_values(by='Authority', ascending=False)[0:30])

	print(Result_df[['PageRank']].sort_values(by='PageRank', ascending=False)[0:30])

	print(Result_df[['SimRank']])
def SortedDict(Dict):
	Dictsorted = sorted(Dict.items())
	sortedlist =list()
	for i in range(len(Dictsorted)):
		sortedlist.append(Dictsorted[i][1])
	#print(sortedlist)
	return sortedlist

def Display_nx_Result(hub,authority,PR):
	
	Result_dict = {
	                'Hub': SortedDict(hub),
	                'Authority': SortedDict(authority),
	                'PageRank' :SortedDict(PR),
	}

	Result_df = pd.DataFrame(Result_dict)
	
	#印出全部矩陣
	pd.set_option('display.max_columns',1000)
	pd.set_option('display.max_rows', 1000) 
	pd.set_option('precision', 3)

	Result_df.sort_values(by='Authority')
	print(Result_df[['Hub','Authority']].sort_values(by='Authority', ascending=False)[0:30])
	print(Result_df[['PageRank']].sort_values(by='PageRank', ascending=False)[0:30])

if __name__ == '__main__':

	# read the input file

	inputFile = 'graph_5.txt'
	edges = edgelist(inputFile)
	m = Nodenum(edges)
	
	A = adjMetrix(edges)
	

	#HITS
	ts = time.clock()
	hub , authority =HITS(A)
	te = time.clock()
	HITStime = str(te-ts)

	#PageRank
	d= 0.15
	ts = time.clock()
	PR = PageRank(A,d)
	te = time.clock()
	PageRanktime = str(te-ts)
	#print(PR)
	
	#SimRank
	ts = time.clock()
	SR =SimRank(A, m)
	te = time.clock()
	SimRanktime = str(te-ts)
	#print(str(SR))

	print('\nGraph Dataset :'+str(inputFile))
	print('\nAdjacency Matrix :\n'+ str(A) )
	#print(hub , authority)
	print('\nImplement >')
	DisplayResult(hub,authority,PR,SR)

	print('\nHITS Algorithm Execution time :'+ HITStime +' sec')
	print('PageRank Algorithm Execution time :'+ PageRanktime+' sec')
	print('SimRank Algorithm Execution time :'+ SimRanktime+' sec')

	G = nx_Graph(edges)
	#DisplayGraph(G)
	G.add_edge(1, 3)
	#G.add_edge(3, 0)

	#DisplayGraph(G)
	#HITS
	nx_ts = time.clock()
	nx_hub , nx_authority =nx_HITS(G)
	nx_te = time.clock()
	nx_HITStime = str(nx_te-nx_ts)
	#print(nx_hub , nx_authority)
	d= 0.15
	#PageRank
	nx_ts = time.clock()
	nx_PR = nx_PageRank(G,1-d)
	nx_te = time.clock()
	nx_PageRanktime = str(nx_te-nx_ts)
	#print(nx_PR)

	print('\nUsing NetworkX >')
	
	Display_nx_Result(nx_hub,nx_authority,nx_PR)

	print('\nHITS Algorithm Execution time :'+ nx_HITStime +' sec')
	print('PageRank Algorithm Execution time :'+ nx_PageRanktime+' sec')
