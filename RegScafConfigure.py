#! /usr/bin/python2

# Copyright © 2021, Mengtian Li, Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import os
import sys
import argparse

descript = "RegScafConfigure.py is used for generate running scripts for RegScaf pipeline.\nConfigure can be run using the command: python RegScafConfigure.py config.txt"
parser = argparse.ArgumentParser(description=descript)

if not os.path.exists("config"):
    os.system("mkdir config")

########## read configuration ##########
with open(sys.argv[1]) as fin_config:
    section = ""
    for line in fin_config:
        line = line.strip()
        if line:
            if line[0] == "#": # section begin
                section_arr = line.split(" ")
                if section_arr[1] == "library_name,":
                    lib_dict = dict()
                    section = "library_name"
                    continue
                elif section_arr[1] == "parameters":
                    section = "parameters"
                    para_dict = dict()
                    # default parameters
                    para_dict["thread"] = 4
                    para_dict["short_seme_penal"] = 1
                    para_dict["long_seme_penal"] = 2
                    para_dict["sspace_link"] = "6"
                    para_dict["depth_threshold"] = 0
                    para_dict["BlasrminMatch"] = 10
                    para_dict["blastfilter_scaffold_size"] = 0
                    para_dict["RegScafPath"] = "."
                    para_dict["reg_m"] = 12
                    para_dict["ScafIter"] = 2
                    para_dict["ContigCorrect"] = False
                    continue

            elif section == "library_name":
                lib_arr = [element.strip() for element in line.split(",")]
                lib_dict[lib_arr[0]] = {"type": lib_arr[1], "insert_size": lib_arr[2], "insert_size_deviation": lib_arr[3], "lib_path1": lib_arr[4], "lib_path2": lib_arr[5]}

            elif section == "parameters":
                para_arr = [element.strip() for element in line.split(",")]
                para_dict[para_arr[0]] = para_arr[1]

PE_arr = []
MP_arr = []
longread_arr = []
for lib_name in lib_dict.keys():
	if lib_dict[lib_name]["type"] == "PE":
		print("PE_Lib",lib_name)
		PE_arr.append(lib_name)
	elif lib_dict[lib_name]["type"] == "MP":
		print("MP_Lib",lib_name)
		MP_arr.append(lib_name)
	elif lib_dict[lib_name]["type"] == "3GS":
		print("LongRead_Lib",lib_name)
		longread_arr.append(lib_name)
########################################

############ BlastFilter.sh ############
with open("config/BlastFilter.sh", "w+") as fout_blastfilter:
    print("dir="+str(para_dict["RegScafPath"])+"/scripts","\n",file=fout_blastfilter)
    if int(para_dict["blastfilter_scaffold_size"]) != 0:
        print("$dir/out_scaffold_long current_reference.fasta length.txt 0",file=fout_blastfilter)
        print("$dir/FilterScaffbyLen_L current_reference.fasta ref_L.fa", para_dict["blastfilter_scaffold_size"],file=fout_blastfilter)
        print("$dir/FilterScaffbyLen_S current_reference.fasta ref_S.fa", para_dict["blastfilter_scaffold_size"],file=fout_blastfilter)
        print("makeblastdb -dbtype nucl -out ref_L_db -in ref_L.fa",file=fout_blastfilter)
        print("blastn -query ref_S.fa -out ref_blastn -db ref_L_db -outfmt 6 -evalue 1e-7 -num_threads ",para_dict["thread"]," -perc_identity 95 -word_size 50",file=fout_blastfilter)
        print("perl $dir/filter_fasta_by_map_rate.pl current_reference.fasta ref_blastn length.txt 0.95 filtered_ref.fa",file=fout_blastfilter)
    else:
        print("ln -s current_reference.fasta filtered_ref.fa",file=fout_blastfilter)
########################################

############ libraries.txt #############
with open("config/libraries.txt", "w+") as fout_libraries:    
    for lib_name in MP_arr:
        if lib_dict[lib_name]["type"] == "PE":
            orient = "FR"
        elif lib_dict[lib_name]["type"] == "MP":
            orient = "RF"
        else:
            print("wrong pair_orientation in", lib_name)
        print("lib"+lib_name, "TAB", "tab"+lib_name+".txt", lib_dict[lib_name]["insert_size"], float(lib_dict[lib_name]["insert_size_deviation"])/int(lib_dict[lib_name]["insert_size"]), orient,file=fout_libraries)
        
########################################

############ RegLib.txt #############
with open("config/RegLib.txt", "w+") as fout_RegLib:
    max_in = 1000
    max_sd = 10
    for lib_name in PE_arr:
        print("lib"+lib_name, lib_dict[lib_name]['type'], "tab"+lib_name+".txt", lib_dict[lib_name]["insert_size"], float(lib_dict[lib_name]["insert_size_deviation"]),file=fout_RegLib)
    for lib_name in MP_arr:
        print("lib"+lib_name, "MP", "tab"+lib_name+".txt", lib_dict[lib_name]["insert_size"], float(lib_dict[lib_name]["insert_size_deviation"]),file=fout_RegLib)
        max_in = max(int(lib_dict[lib_name]["insert_size"]),max_in)
        max_sd = max(int(lib_dict[lib_name]["insert_size_deviation"]),max_sd)
    for lib_name in longread_arr:
        print("lib"+lib_name, "3GS","3GSTab"+lib_name, 0, 0,file=fout_RegLib)
########################################

############## AllReg-iter.sh ###############
with open("config/AllReg-iter.sh", "w+") as fout_allreg:
    print("#!/bin/bash\n",file=fout_allreg)
    print("#directory to scripts",file=fout_allreg)
    print("dir="+str(para_dict["RegScafPath"])+"/scripts","\n",file=fout_allreg)
    print("python $dir/PreFilterTab.py detail.txt RegLib.txt \npython $dir/ReCalLibInsertSizeDF.py detail.txt RegLib.txt\n",file=fout_allreg)
    print("python $dir/ContigGraph-LTS-multiMean.py -d detail.txt -l RegLib.txt -m",para_dict["reg_m"],"-s 300 -M 100 -t",para_dict["thread"],"> ConnectGraph.log 2>err\n",file=fout_allreg)
    print("iter=1\nwhile ((iter<"+str(para_dict["ScafIter"])+"))\ndo\n\tif [ -d \"./final_dir\" ]; then\n\t\tmv final_dir Final_dir\n\tfi\n\tpython $dir/SplitFinal_all.py -d detail.txt -sd $dir/ -m",int(int(para_dict["reg_m"])/2),"-mi",max_in,"-ms",max_sd,"-t",para_dict["thread"]," > SplitFinal.log_iter$iter",file=fout_allreg)
    print("\tpython $dir/ScaffoldGraph-LTS-multimean.py -d detail.txt -l RegLib.txt -m",para_dict["reg_m"],"-M $((100+$iter*50)) -t",para_dict["thread"]," > ScafGraph_iter$iter.log 2>>err",file=fout_allreg)
    print("\tmv Final_dir Final_iter$iter\n\tmv SplitFinal SplitFinal_iter$iter\n\t((iter++))\ndone\n",file=fout_allreg)
    if para_dict["ContigCorrect"]:
        print("python $dir/MergeFinal.py -d detail.txt -sd $dir/ -t",para_dict["thread"],"-m",int(int(para_dict["reg_m"])/2),"-ms",max_sd,"-mi",max_in,"-o regScaf -C >Merge.log 2>>err \npython $dir/MakeEvidence-noSingleG.py evidence.reg\n",file=fout_allreg)
    else:
        print("python $dir/MergeFinal.py -d detail.txt -sd $dir/ -t",para_dict["thread"],"-m",int(int(para_dict["reg_m"])/2),"-ms",max_sd,"-mi",max_in,"-o regScaf >Merge.log 2>>err \npython $dir/MakeEvidence-noSingleG.py evidence.reg\n",file=fout_allreg)
    print("$dir/GenomeStat regScaf_scaf.fa",file=fout_allreg)
########################################

'''
############## AllReg.sh ###############
with open("config/AllReg.sh", "w+") as fout_allreg:
    print("#!/bin/bash\n",file=fout_allreg)
    print("#directory to scripts",file=fout_allreg)
    print("dir="+str(para_dict["RegScafPath"]),"\n",file=fout_allreg)
    print("python $dir/PreFilterTab.py detail.txt RegLib.txt \npython $dir/ReCalLibInsertSizeDF.py detail.txt RegLib.txt\n",file=fout_allreg)
    print("python $dir/ContigGraph-LTS-multiMean.py -d detail.txt -l RegLib.txt -m",para_dict["reg_m"],"-s 300 -M 100 -t",para_dict["thread"],"> ConnectGraph.log 2>err \npython $dir/SplitFinal_all.py -d detail.txt -sd $dir/ -m",int(int(para_dict["reg_m"])/2),"-mi",max_in,"-t",para_dict["thread"]," > SplitFinal.log\n",file=fout_allreg)
    print("python $dir/ScaffoldGraph-LTS-multimean.py -d detail.txt -l RegLib.txt -m",para_dict["reg_m"],"-M 150 -t",para_dict["thread"]," > ScafGraph.log 2>>err\n",file=fout_allreg)
    print("python $dir/MergeFinal.py -d detail.txt -sd $dir/ -t",para_dict["thread"],"-m",int(int(para_dict["reg_m"])/2),"-ms",max_sd,"-mi",max_in,"-o regScaf >Merge.log 2>>err \npython $dir/MakeEvidence-noSingleG.py evidence.reg",file=fout_allreg)
########################################
'''

############## SSPACE.sh ###############
with open("config/SSPACE.sh", "w+") as fout_sspace:
    print("SSPACE_Standard_v3.0.perc10perc90.pl -l libraries.txt -s initial_contig.fa -x 0 -k", para_dict["sspace_link"], "-a 0.7 -n 15 -T", para_dict["thread"]," -S 0 -b sspace_k"+para_dict["sspace_link"], " > /dev/null 2>&1",file=fout_sspace)
    print("cp sspace_k"+para_dict["sspace_link"]+"/sspace_k"+para_dict["sspace_link"]+".final.evidence ./evidence.txt.sspace",file=fout_sspace)
    print("cp ~/BAUM/quantileCode/ParseSSPACEEvidence.perc10perc90.pl .",file=fout_sspace)
    print("perl ParseSSPACEEvidence.perc10perc90.pl evidence.txt.sspace evidence.txt",file=fout_sspace)
########################################


############## MakeTab.sh ###############
with open("config/MakeTab.sh", "w+") as fout_tab:
    print("dir="+str(para_dict["RegScafPath"])+"/scripts","\n",file=fout_tab)
    for lib_name in longread_arr: ##Filter hanging out
        print("awk '($7<200 || $10<200) && (($9-$8)<200 ||($12-$11)<200)' BLASROut"+lib_name+"> FBLASROut"+lib_name,file=fout_tab)
    for lib_name in PE_arr:
        print("$dir/MakeTabFile_v2.0", lib_name+"_R1.sam", lib_name+"_R2.sam", "initial_contig.fa No No highdepth.bed", "tab"+lib_name+".txt", "&",file=fout_tab)
    for lib_name in MP_arr:
        print("$dir/MakeTabFile_v2.0", lib_name+"_R1.sam", lib_name+"_R2.sam", "initial_contig.fa No No highdepth.bed", "tab"+lib_name+".txt", "&",file=fout_tab)
    for lib_name in longread_arr:
        print("python $dir/BlasrIntoLink.py FBLASROut"+lib_name+" 3GSTab"+lib_name, "&",file=fout_tab)
    print("wait",file=fout_tab)

########################################

############### Mapping.sh ################
with open("config/Mapping.sh", "w+") as fout_map:
    print("dir="+str(para_dict["RegScafPath"])+"/scripts","\n",file=fout_map)
    print("if [ -d \"index\" ]; then \n  mv index index_old \nfi\n$dir/sort103 reference.fasta index",file=fout_map)
    for lib_name in PE_arr:
        print("$dir/seme_v2.0 index/", lib_dict[lib_name]["lib_path1"], "-o", lib_name+"_R1.sam", "-c 20 -s 15 -ml 50 -w 2 -penal", para_dict["short_seme_penal"], "-t", para_dict["thread"], "-rc seme_mapresult_"+lib_name+"_R1.txt",file=fout_map)
        print("$dir/seme_v2.0 index/", lib_dict[lib_name]["lib_path2"], "-o", lib_name+"_R2.sam", "-c 20 -s 15 -ml 50 -w 2 -penal", para_dict["short_seme_penal"], "-t", para_dict["thread"], "-rc seme_mapresult_"+lib_name+"_R2.txt",file=fout_map)
    for lib_name in MP_arr:
        print("$dir/seme_v2.0 index/", lib_dict[lib_name]["lib_path1"], "-o", lib_name+"_R1.sam", "-c 20 -s 15 -ml 50 -w 2 -penal", para_dict["long_seme_penal"], "-t", para_dict["thread"], "-rc seme_mapresult_"+lib_name+"_R1.txt",file=fout_map)
        print("$dir/seme_v2.0 index/", lib_dict[lib_name]["lib_path2"], "-o", lib_name+"_R2.sam", "-c 20 -s 15 -ml 50 -w 2 -penal", para_dict["long_seme_penal"], "-t", para_dict["thread"], "-rc seme_mapresult_"+lib_name+"_R2.txt",file=fout_map)
    for lib_name in longread_arr:
        print("blasr", lib_dict[lib_name]["lib_path1"], "initial_contig.fa --minMatch",para_dict['BlasrminMatch'],"--minPctSimilarity 85 --maxScore -600 --minReadLength 1000 --minSubreadLength 1000 --minAlnLength 100 --nproc",para_dict["thread"],"--out BLASROut"+lib_name,file=fout_map)
########################################


############### Bed.sh #################
with open("config/Bed.sh", "w+") as fout_bedsh:
    print("dir="+str(para_dict["RegScafPath"])+"/scripts","\n",file=fout_bedsh)
    print("mode=`$dir/CalPrimDepth reference.fasta $1 PrimDepth.txt $2| tail -1 | awk '{print($2*1.8)}'`",file=fout_bedsh)
    if para_dict["depth_threshold"] == 0:
        print("$dir/HighDepthRegion PrimDepth.txt highdepth.bed $mode",file=fout_bedsh)
    else:
        print("$dir/HighDepthRegion PrimDepth.txt highdepth.bed", para_dict["depth_threshold"],file=fout_bedsh)
    print("wait",file=fout_bedsh)
########################################


############ name_list.txt #############
with open("config/name_list.txt", "w+") as fout_namelist:
    for lib_name in PE_arr:
        print(lib_name+"_R1.sam\n", lib_name+"_R2.sam",file=fout_namelist)
    for lib_name in MP_arr:
        print(lib_name+"_R1.sam\n", lib_name+"_R2.sam",file=fout_namelist)
########################################
