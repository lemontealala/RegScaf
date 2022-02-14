from sklearn.neighbors import KernelDensity
from scipy import signal
import networkx as nx
import numpy as np

def ClusterPeak(List,sd):
	sd = max(10,sd)
	x=np.array(sorted(List))
	y=[[i] for i in x]
	kde = KernelDensity(kernel='gaussian', bandwidth=min(sd/2,100)).fit(y)
	log_density = kde.score_samples(y)
	est_density = np.exp(log_density)
	num_peak_3 = signal.find_peaks(est_density,distance=3)
	peaks = [x[i] for i in num_peak_3[0]]
	#print("Peaks:",peaks)
	peak_density = [est_density[i] for i in num_peak_3[0]]
	border = [(peaks[i]+peaks[i+1])/2 for i in range(len(peaks)-1)]
	##Those jumping point should also be border
	i = 0
	while i<len(x)-1:
		if x[i+1]-x[i]>sd :
			border.append(x[i]+sd/2)
		i +=1
	clusters = []
	ai = min(x)-1
	border.sort()
	#print("Borders:",border)
	i = 0
	while i < len(border):
		bi = border[i]
		group = list(x[(ai<x)*(x<=bi)])
		if len(group)>0:
			clusters.append(group)
		ai = bi
		i+=1
	clusters.append(list(x[ai<x]))
	if len(clusters)<2:
		#print('Total:',len(List),'the number of cluster is ' + str(len(clusters)))
		return clusters
	else:
		i = 0
		while i<len(clusters)-1:
			#print("Diff between clusters:",max(clusters[i]),min(clusters[i+1]))
			if min(clusters[i+1])-max(clusters[i])<sd/2:#100:
				#print("Diff between clusters <sd/2 ! Merge them.",len(clusters))
				clusters[i]+=clusters[i+1]
				del clusters[i+1]
				#print(len(clusters))
			else:
				#print(clusters[i],clusters[i+1])
				i+=1
	#print('Total:',len(List),'the number of cluster is ' + str(len(clusters)))
	return clusters

		
def Select_Col(select_row,Coff_mat,Weight_mat,current_ctg,truc):
	select_col = set()
	G=nx.Graph()
	for i in current_ctg:
		G.add_node(contig_list[i])
	pair = {}
	for row in select_row[1:]:
		line=list(Coff_mat[row,:])
		select_col.add(line.index(1))
		select_col.add(line.index(-1))
		n1=current_ctg[line.index(1)]
		n2=current_ctg[line.index(-1)]
		G.add_edge(contig_list[n1],contig_list[n2])
	allSubG = nx.connected_component_subgraphs(G)
	IsConnect = 0
	for subG in allSubG:
		print("Here subG.size",len(subG.nodes()))
		if contig_list[current_ctg[0]] in subG.nodes() and len(subG.nodes())==len(current_ctg):		
			IsConnect = 1
			print("Here IsConnect!")
		else:
			subG_nodes = list(set(subG.nodes()))
			ProRegFunc(subG_nodes)
	return(IsConnect,list(select_col))
	
def RegProc(Coff_mat,Weight_mat,y,Cov_mat):
	doubleCoff_0 = np.dot(Coff_mat.T,Weight_mat)
	doubleCoff = np.dot(doubleCoff_0,Coff_mat)
	Solution_mat_0 = np.dot(np.linalg.pinv(doubleCoff),Coff_mat.T)
	Lamba_mat = np.dot(Solution_mat_0,Weight_mat)
	Solution_mat = np.dot(Lamba_mat, y.T)
	y_predict = np.dot(Coff_mat,Solution_mat)
	y_diff = y_predict - y
	abs_error = abs(y_diff)
	Solution_cov = np.dot(np.dot(Lamba_mat,Cov_mat),Lamba_mat.T)
	return (Solution_mat,abs_error,Solution_cov)

def RegProcSol(Coff_mat,Weight_mat,y): #,Cov_mat
	doubleCoff_0 = np.dot(Coff_mat.T,Weight_mat)
	doubleCoff = np.dot(doubleCoff_0,Coff_mat)
	Solution_mat_0 = np.dot(np.linalg.pinv(doubleCoff),Coff_mat.T)
	Lamba_mat = np.dot(Solution_mat_0,Weight_mat)
	Solution_mat = np.dot(Lamba_mat, y.T)
	#Solution_cov = np.dot(np.dot(Lamba_mat,Cov_mat),Lamba_mat.T)
	return (Solution_mat) #,Solution_cov

def SolCovFinal(Coff_mat, Weight_mat, Cov_mat):
	doubleCoff_0 = np.dot(Coff_mat.T,Weight_mat)
	doubleCoff = np.dot(doubleCoff_0,Coff_mat)
	Solution_mat_0 = np.dot(np.linalg.pinv(doubleCoff),Coff_mat.T)
	Lamba_mat = np.dot(Solution_mat_0,Weight_mat)
	Solution_cov = np.dot(np.dot(Lamba_mat,Cov_mat),Lamba_mat.T)
	return Solution_cov

def findDemarcation(r_list,maxerror):
    size = len(r_list)
    if size<3:
        return min(r_list)
    s_list = sorted(r_list)
    mid_id = int(size*3/4)
    k = mid_id+1
    while k<size-1:
        last = s_list[k]
        diff = s_list[k+1]-last
        if last<maxerror: 
            k +=1
            continue
        if diff/last>3 and s_list[k+1]>1000: #r>10*maxerror
            print("mid",mid_id,"Found Demarcation:",size,k,last,s_list[k+1])
            return(last)
        else:
            k +=1
    if size>200:
        return (np.percentile(r_list,95)-1)
    else:
        return (np.percentile(r_list,98)-1)

def FindDemarcation(r_list,maxerror):
    size = len(r_list)
    if size<3:
        return min(r_list)
    s_list = sorted(r_list)
    mid_id = int(size*3/4)+1
    last = s_list[mid_id-1]
    for r in s_list[mid_id:]:
        diff = r-last
        if last==0: 
            last = r
            continue
        if diff/last>3 and r>1000: #r>maxerror
            print("Found Demarcation:",size,s_list.index(last),last,r)
            '''if size>1000:
                print("Here1000:")
                np.savetxt("FindDema_example_1000.txt",s_list,fmt="%d")
            elif size>500:
                print("Here500:")
                np.savetxt("FindDema_example_500.txt",s_list,fmt="%d")
            elif size>100:
                print("Here100:")
                np.savetxt("FindDema_example_100.txt",s_list,fmt="%d")'''
            return(last)
        last = r
    if size>200:
        return np.percentile(r_list,95)#95)
    else:
        return np.percentile(r_list,98)
	
def TrimMean(x):
	lower_bound = np.percentile(x,5)
	upper_bound = np.percentile(x,95)
	r_list = [i for i in x if (i>=lower_bound and i<=upper_bound)]	#list(filter(lambda y:y>=lower_bound and y<=upper_bound, list))
	if len(r_list)>1:
		return np.median(r_list) #np.mean
	else:
		return np.median(x)
	
def MADN(list):
	list=np.array(list)
	MAD = np.median(abs(list-np.median(list)))
	return MAD/0.675
	
def MedianStd(list):
	list=np.array(list)
	size=len(list)
	if size<2: return 0
	#mean=np.mean(list1)
	median=np.median(list)
	m_list=median*np.ones(size)
	Std=sum((list-m_list)**2)/size
	return sqrt(Std/2)
	
def del_col_and_update(Coff_mat,Weight_mat,y,Cov_mat,subG,current_ctg):
	global contig_list
	G_col = []
	new_current = []
	for node in subG.nodes():
		id_in_contigl = contig_list.index(node)
		id_in_current = current_ctg.index(id_in_contigl)
		G_col.append(id_in_current)
		new_current.append(id_in_contigl)
	Coff_mat = Coff_mat[:,G_col]
	j=1
	while j <= Coff_mat.shape[0]-1:
		line=list(Coff_mat[j,:])
		if (1 not in line) or (-1 not in line):	#row line should be deleted
			Coff_mat=np.delete(Coff_mat,j,axis=0)
			y=np.delete(y,j)
			Cov_mat=np.delete(Cov_mat,j,axis=0)
			Cov_mat=np.delete(Cov_mat,j,axis=1)
			Weight_mat=np.delete(Weight_mat,j,axis=0)
			Weight_mat=np.delete(Weight_mat,j,axis=1)
		else : j +=1
	return(Coff_mat,Cov_mat,Weight_mat,new_current,y)	