import re
from Bio import SeqIO
import numpy as np
import os

def writeResults(f, number, length, contig):
    f.write(">contig%d\n"%number)
    line=int(length/100)
    for k in range(line):
        f.write(contig[100*k:100*(k+1)]+'\n')
    if length>line*100:
        f.write(contig[100*line:]+'\n')

fout = open('initial_contig.fa','w')
fout_length = open('detail.txt','w')
fout_evidence = open('evidence_ini.txt','w')

if (os.path.exists("CONTIG")):
	setdir_com="rm -r CONTIG"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir CONTIG"
print(setdir_com)
os.system(setdir_com)

contigNum = 0
scafNum = 1
for seq_record in SeqIO.parse('genome.fasta', 'fasta'):
    scaf_content=str(seq_record.seq)    
    scaf_break_arr = re.split('(N+)', scaf_content)
    contig_count = int(np.ceil(len(scaf_break_arr)/2))
    fout_evidence.write("~~~~~\n%s\t%d\t%d\n"%(scafNum,len(scaf_content), contig_count))
    scafNum += 1
    for i in range(contig_count):
        new_contig = scaf_break_arr[2*i]
        if i==contig_count-1:
            scaf_N_count = 0
        else:
            scaf_N_count = len(scaf_break_arr[2*i+1])
        length=len(new_contig)
        contigNum+=1
        fout_length.write('contig%d\t1\t%d\n'%(contigNum, length))
        fout_evidence.write('%d\tf\t%d\t0\t%d\t0\n'%(contigNum, length, scaf_N_count))
        writeResults(fout, contigNum, length, new_contig)
        writeResults(open("CONTIG/contig"+str(contigNum)+".fasta",'w'), contigNum, length, new_contig)

fout.close()
fout_length.close()
fout_evidence.close()
