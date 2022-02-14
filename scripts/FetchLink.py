
def GetLongLink(BLASR_line):
	pairline=BLASR_line.strip().split(" ")
	c1 = int(pairline[0])
	c2 = int(pairline[1])
	d = -2*(int(pairline[2]))+1 #0-->1,1-->-1
	diff = int(pairline[3])
	return(c1,c2,diff,d)

def CalMPDiff(tabline,InsertSize):# 	——> ingap <——
	pairline = tabline.strip().split("\t")
	c1=int(pairline[0][6:])
	c2=int(pairline[3][6:])
	p1 = int(pairline[1])
	m1 = int(pairline[2])
	p2 = int(pairline[4])
	m2 = int(pairline[5])
	ingap = InsertSize-abs(p1-m1)-abs(p2-m2)
	if p1>m1:
		d1 = 1	#consistent with paired-end read
		x1 = -p1
	elif p1<m1:
		d1 =0	#inconsistent
		x1 = p1
	if p2<m2:
		d2 = 1	#consistent with paired-end read
		x2 = ingap-p2
	elif p2>m2:
		d2 =0	#inconsistent
		x2 = ingap+p2
	diff = x2-x1
	if d1==0:
		diff = -diff
	if d1==d2:
		return(c1,c2,diff,1)
	else:
		return(c1,c2,diff,-1)

def CalPEDiff(tabline,InsertSize): # <—— ingap ——> sign of difference consistent with the prior contig
	pairline = tabline.strip().split("\t")
	c1=int(pairline[0][6:])
	c2=int(pairline[3][6:])
	p1 = int(pairline[1])
	m1 = int(pairline[2])
	p2 = int(pairline[4])
	m2 = int(pairline[5])
	ingap = InsertSize#+abs(p1-m1)+abs(p2-m2)
	if p1<m1:
		d1 = 1	#consistent with paired-end read
		x1 = -p1
	elif p1>m1:
		d1 =0	#inconsistent
		x1 = p1
	if p2>m2:
		d2 = 1	
		x2 = ingap-p2
	elif p2<m2:
		d2 =0	
		x2 = ingap+p2
	diff = x2-x1
	if d1==0:	#sign of difference must be in accordance with the prior contig's direction
		diff = -diff
	if d1==d2:
		return(c1,c2,diff,1)
	else:
		return(c1,c2,diff,-1)
