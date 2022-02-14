#SortScaffold.py
import sys

ScafFile = open(sys.argv[1],'r')
ScafFile_new = open(sys.argv[2],'w')
scaf_num = 1
for line in ScafFile:
	if line=='\n':
		continue
	if line[0]==">":
		line_new = ">scaffold"+str(scaf_num)
		print(line_new,file = ScafFile_new)
		scaf_num+=1
	else:
		print(line,end="",file = ScafFile_new)
	
ScafFile.close()
ScafFile_new.close()