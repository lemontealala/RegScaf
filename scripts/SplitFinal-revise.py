#python3 MergeAndWriteNewScaffold-revise.py Final.txt truc max_insert max_sd 
import sys
import os
import numpy as np
import copy

def IfOverlap(s1,e1,s2,e2):
	if (s1<s2 and e1-s2>3*max_sd) or (s2<s1 and e2-s1>3*max_sd):
		return True
	else:
		return False

def IfConflict(single_ctg,Scaf):
	for s_ctg in Scaf:
		if IfOverlap(s_ctg['start'],s_ctg['end'],single_ctg['start'],single_ctg['end']):
			return True
	return False

def GetLink(c1,c2):
	global linkdict
	link = 0
	k = str(c1)+"-"+str(c2)
	if k in linkdict:
		link = linkdict[k]
	return link

def GetMostLink(c1,ctg_dict):
	max_link = 0
	max_c2 = c1
	for ctg in ctg_dict:
		link = GetLink(c1,ctg['contig'])
		if link>max_link:
			max_link = link
			max_c2 = ctg['contig']
	return(max_c2,max_link)

def ReCalGap(ctg_dict):
	ctg_dict.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	ctg_num = len(ctg_dict)
	i = 0
	contig_dict = copy.deepcopy(ctg_dict)
	while i < ctg_num-1:
		gap_new = contig_dict[i+1]['start']-contig_dict[i]['end']
		contig_dict[i]['gap'] = int(gap_new)
		i +=1
	contig_dict[i]['gap'] = 0
	return contig_dict
	
def WriteSplitFinal(FinalScaf):
	FinalScaf.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	size = len(FinalScaf)
	scafmaxlen = FinalScaf[size-1]['end']-FinalScaf[0]['start']
	fout_Final = open("SplitFinal/Final_"+str(FinalScaf[0]['contig'])+".txt",'w+')
	print(("longest position is %.1f" % scafmaxlen),file = fout_Final)
	for k in range(len(FinalScaf)):
		dict=FinalScaf[k]
		if k== len(FinalScaf)-1:
			gap=0
		else:	
			next=FinalScaf[k+1]
			gap=next['start']-dict['end']
		print(dict['contig'], '%.1f' % dict['start'], '%.1f' % dict['end'],	dict['direct'],'%.1f' % dict['std'],'%.1f' % gap,'%d' % dict['len'], file=fout_Final)
	fout_Final.close()

##global variable
max_insert = int(sys.argv[3])
max_sd = int(sys.argv[4])

##Get Information 
result = os.popen('tail -1 detail.txt')
res = result.read().split()
ContigCount = int(res[0][6:])

linkFile = open("link_dict",'r')
a = linkFile.read()
linkdict = eval(a)
linkFile.close()	
RegOut = open(sys.argv[1],'r')
truc = int(sys.argv[2])
Reglines = RegOut.readlines()
contig_dict=[]
for line in Reglines[1:]:
	record = line.strip().split(' ')
	contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6])),'newctg':0})
ctg_num = len(contig_dict)	

##Split into conflict-free strands
contig_dict.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
AllCtg_dict = copy.deepcopy(contig_dict)
Split_Scaf = []
Split_Scaf_index = []
selected_index = []
scaf_num = 0
CtgScaf_dict = {}
remained_index = [ x for x in range(ctg_num)]
while len(remained_index)>0:
	Select_Scaf = []
	i = remained_index[0]
	selected_index = [remained_index[0]]
	CtgScaf_dict[AllCtg_dict[i]['contig']] = scaf_num
	Select_Scaf.append(AllCtg_dict[i])	
	k_i = 0
	k_j = 0
	while k_i<len(remained_index)-1 :
		#print("Now i:",AllCtg_dict[i]['contig'])
		window_start = AllCtg_dict[i]['end']-3*max_sd
		window_end = window_start+max_insert+6*max_sd
		max_link = 0
		max_j = i
		k_j = k_i+1
		j = remained_index[k_j]
		while k_j<len(remained_index) and AllCtg_dict[j]['start'] <window_end:#AllCtg_dict[j]['start'] >window_start and
			j = remained_index[k_j]
			c1 = AllCtg_dict[i]['contig']
			c2 = AllCtg_dict[j]['contig']
			link = GetLink(c1,c2)
			if link>max_link:
				max_link = link
				max_j = j
			k_j +=1
		if max_link>truc:
			print("max_link",max_link,AllCtg_dict[i]['contig'],AllCtg_dict[max_j]['contig'])
			selected_index.append(max_j)
			CtgScaf_dict[AllCtg_dict[max_j]['contig']] = scaf_num
			Select_Scaf.append(AllCtg_dict[max_j])
			i = max_j
			k_i = remained_index.index(i)
		else:
			break
	Split_Scaf.append(Select_Scaf)
	Split_Scaf_index.append(selected_index)
	scaf_num +=1
	remained_index = sorted(list(set(remained_index)-set(selected_index)))
	
print("Before join:",[[scaf_ctg['contig'] for scaf_ctg in Scaf] for Scaf in Split_Scaf])
##Join single scaffold into their most linked set if no conflict
k = 0
rm_id = []
while k<len(Split_Scaf):
	Select_Scaf = Split_Scaf[k]
	if len(Select_Scaf)<2:
		(max_c2,max_l) = GetMostLink(Select_Scaf[0]['contig'],AllCtg_dict)
		#print("Single Scaf:",Select_Scaf[0]['contig'],max_c2,max_l)
		if max_l>truc :
			id_l = CtgScaf_dict[max_c2]
			if IfConflict(Select_Scaf[0],Split_Scaf[id_l]):
				print("Conflict if join!",max_c2,max_l,Select_Scaf[0]['contig'],id_l)
				k+=1
				continue
			print("Join:",max_c2,max_l,Select_Scaf[0]['contig'],id_l)
			Split_Scaf[id_l] += Select_Scaf
			for ctg in Select_Scaf: #Update the dictionary
				CtgScaf_dict[ctg['contig']] = id_l
			rm_id.append(k)
			#Split_Scaf.remove(Select_Scaf)
			#k -=1
	k +=1
remain_id = list(set(range(len(Split_Scaf)))-set(rm_id))
FSplit_Scaf = [Split_Scaf[i] for i in remain_id]
print("Final:",[[scaf_ctg['contig'] for scaf_ctg in Scaf] for Scaf in FSplit_Scaf])

for Scaf in FSplit_Scaf:
	WriteSplitFinal(Scaf)