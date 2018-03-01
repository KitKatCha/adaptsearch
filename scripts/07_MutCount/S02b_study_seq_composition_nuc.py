#!/usr/bin/env python
## Author: Eric FONTANILLAS
## Date: 21.12.10
## Last Version : 12/2017 by Victor Mataigne
## Object: Test for compositional bias in genome and proteome as marker of thermal adaptation (comparison between 2 "hot" species: Ap and Ps and one "cold" species: Pg)

import sys,os,shutil,subprocess, string, itertools
from functions import simplify_fasta_name, dico

##################
###### DEF2 ######
##################
def base_composition(seq):
    count_A=string.count(seq, "A") + string.count(seq, "a")
    count_T=string.count(seq, "T") + string.count(seq, "t")
    count_C=string.count(seq, "C") + string.count(seq, "c")
    count_G=string.count(seq, "G") + string.count(seq, "g")
    ## 3 ## Nucleotide proportions
    ln = count_C+count_G+count_T+count_A
    if (ln!=0):

        CG = count_C+count_G
        AT = count_T+count_A
      
        AG = count_A+count_G
        TC = count_T+count_C

        ## 1 ## Search for compositional bias in genome as marker of thermal adaptation: CG vs AT
        ratio_CG_AT = float(CG)/float(AT)
        percent_CG = float(CG)/(float(AT) + float(CG))*100
  
        ## 2 ## Search for compositional bias in genome as marker of thermal adaptation: AG vs TC
        ratio_purine_pyrimidine=float(AG)/float(TC)
        percent_purine=float(AG)/(float(AG)+float(TC))*100      
  
 
        prop_A = float(count_A)/float(ln)
        prop_T = float(count_T)/float(ln)
        prop_C = float(count_C)/float(ln)
        prop_G = float(count_G)/float(ln)
    else:
        percent_CG=0
        percent_purine=0
        prop_A=0
        prop_T=0
        prop_C=0
        prop_G=0
      
    return(percent_CG, percent_purine, prop_A, prop_T, prop_C, prop_G)
##############################################

##################
###### DEF3 ######
##################
def purine_loading(seq):
    count_A=string.count(seq, "A") + string.count(seq, "a")
    count_T=string.count(seq, "T") + string.count(seq, "t")
    count_C=string.count(seq, "C") + string.count(seq, "c")
    count_G=string.count(seq, "G") + string.count(seq, "g")
    ## 3 ## Nucleotide proportions
    TOTAL = count_C+count_G+count_T+count_A
    if (TOTAL!=0):
     
        ## PLI : Purine loading indice (Forsdyke et al.)
        # (G-C)/N * 1000 et (A-T)/N * 1000
  
        DIFF_GC = count_G - count_C
        DIFF_AT = count_A - count_T

        # Per bp
        PLI_GC = float(DIFF_GC)/float(TOTAL)
        PLI_AT = float(DIFF_AT)/float(TOTAL)
  
        # Per 1000 bp
        PLI_GC_1000 = PLI_GC*1000
        PLI_AT_1000 = PLI_AT*1000
    else:
        DIFF_GC=0
        DIFF_AT=0
        PLI_GC=0
        PLI_AT=0
        PLI_GC_1000=0
        PLI_AT_1000=0

    return(TOTAL, DIFF_GC, DIFF_AT,PLI_GC,PLI_AT,PLI_GC_1000,PLI_AT_1000)
##############################################

###################
### RUN RUN RUN ###
###################

##Create specific folders
Path_IN_loci_NUC = "./IN_NUC"
outpath= "./OUT"
os.makedirs(Path_IN_loci_NUC)
os.makedirs(outpath)

infiles = str.split(sys.argv[1], ",")
for file in infiles:
    if os.stat(file).st_size == 0:
        pass
    else:
        os.system("cp %s %s" %(file, Path_IN_loci_NUC))

## 1 ## List taxa
LT=[]
cmd="grep '>' %s" % sys.argv[2]
result = subprocess.check_output(cmd, shell=True)
result=result.split('\n')
for i in result:   
    sp=i[1:]
    if sp !='':
        LT.append(sp)
print LT

#LT = ["Ap", "Ac", "Pg"]
#os.system("grep '>' %s" %(sys.argv[1]))

## 2 ## PathIN
# fileIN_properties = open("01_AminoAcid_Properties2.csv", "r")
Lloci_NUC = os.listdir(Path_IN_loci_NUC)


## 3 ## PathOUT
## 3.1 ## NUC composition
fileOUT_NUC=open("./OUT/nuc_compositions.csv","w")
fileOUT_NUC.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_NUC.write("%s_prop_A,%s_prop_T,%s_prop_C,%s_prop_G," %(taxa,taxa,taxa,taxa))
fileOUT_NUC.write("%s_prop_A,%s_prop_T,%s_prop_C,%s_prop_G" %(LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_NUC.write("\n")

## 3.2 ## NUC percent_GC
fileOUT_percent_GC=open("./OUT/percent_GC.csv","w")
fileOUT_percent_GC.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_percent_GC.write("%s_percent_GC," %(taxa))
fileOUT_percent_GC.write("%s_percentGC" %(LT[-1]))
fileOUT_percent_GC.write("\n")

## 3.3 ## NUC percent_purine
fileOUT_percent_purine=open("./OUT/percent_purine.csv","w")
fileOUT_percent_purine.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_percent_purine.write("%s_percent_purine," %(taxa))
fileOUT_percent_purine.write("%s_percent_purine" %(LT[-1]))
fileOUT_percent_purine.write("\n")

## 3.4 ## Purine Load
fileOUT_Purine_Load=open("./OUT/Purine_Load_Indice.csv", "w")
fileOUT_Purine_Load.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_Purine_Load.write("%s_TOTAL,%s_DIFF_GC,%s_DIFF_AT,%s_PLI_GC1000,%s_PLI_AT1000," %(taxa,taxa,taxa,taxa,taxa))
fileOUT_Purine_Load.write("%s_TOTAL,%s_DIFF_GC,%s_DIFF_AT,%s_PLI_GC1000,%s_PLI_AT1000" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_Purine_Load.write("\n")

#####################
## 4 ## Process Loci
#####################
for locus in Lloci_NUC:
    if not os.path.isfile("%s/%s" %(Path_IN_loci_NUC, locus)) :
        pass
    else :
        print locus
        path_locus = "%s/%s" %(Path_IN_loci_NUC, locus)
        bash = dico(path_locus,LT) 

        fileOUT_NUC.write("%s," %locus)
        fileOUT_percent_GC.write("%s," %locus)
        fileOUT_percent_purine.write("%s," %locus)
        fileOUT_Purine_Load.write("%s," %locus)

        for taxa in LT[0:-1]:    
            if taxa in bash.keys():        
                seq = bash[taxa]            
                percent_GC, percent_purine,prop_A, prop_T, prop_C, prop_G = base_composition(seq)   ### DEF2 ###
                TOTAL, DIFF_GC, DIFF_AT,PLI_GC,PLI_AT,PLI_GC_1000,PLI_AT_1000 = purine_loading(seq) ### DEF3 ###
                fileOUT_NUC.write("%.5f,%.5f,%.5f,%.5f," %(prop_A,prop_T,prop_C,prop_G))
                fileOUT_percent_GC.write("%.5f," %percent_GC)
                fileOUT_percent_purine.write("%.5f," %percent_purine)
                fileOUT_Purine_Load.write("%d,%d,%d,%.5f,%.5f," %(TOTAL, DIFF_GC, DIFF_AT,PLI_GC_1000, PLI_AT_1000))
            else:
                fileOUT_NUC.write("%s,%s,%s,%s," %("NA","NA","NA","NA"))
                fileOUT_percent_GC.write("%s," %"NA")
                fileOUT_percent_purine.write("%s," %"NA")
                fileOUT_Purine_Load.write("%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA"))
          
        if LT[-1] in bash.keys():
            seq = bash[LT[-1]]            
            percent_GC, percent_purine,prop_A, prop_T, prop_C, prop_G = base_composition(seq)   ### DEF2 ###
            TOTAL, DIFF_GC, DIFF_AT,PLI_GC,PLI_AT,PLI_GC_1000,PLI_AT_1000 = purine_loading(seq) ### DEF3 ###
            fileOUT_NUC.write("%.5f,%.5f,%.5f,%.5f" %(prop_A,prop_T,prop_C,prop_G))
            fileOUT_percent_GC.write("%.5f" %percent_GC)
            fileOUT_percent_purine.write("%.5f" %percent_purine)
            fileOUT_Purine_Load.write("%d,%d,%d,%.5f,%.5f" %(TOTAL, DIFF_GC, DIFF_AT,PLI_GC_1000, PLI_AT_1000))
        else:
            fileOUT_NUC.write("%s,%s,%s,%s" %("NA","NA","NA","NA"))
            fileOUT_percent_GC.write("%s" %"NA")
            fileOUT_percent_purine.write("%s" %"NA")
            fileOUT_Purine_Load.write("%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA"))


        fileOUT_NUC.write("\n")
        fileOUT_percent_GC.write("\n")
        fileOUT_percent_purine.write("\n")
        fileOUT_Purine_Load.write("\n")
fileOUT_NUC.close()
fileOUT_percent_GC.close()
fileOUT_percent_purine.close()
fileOUT_Purine_Load.close()
