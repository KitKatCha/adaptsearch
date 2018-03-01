#!/usr/bin/env python
# -*- coding: ascii -*-
## Author: Eric FONTANILLAS
## Date: 21.12.10
## Last Version : 12/2017 by Victor Mataigne
## Object: Test for compositional bias in genome and proteome as marker of thermal adaptation (comparison between 2 "hot" species: Ap and Ps and two "cold" species: Pg, Pp)

import sys,os,shutil,subprocess,string, itertools
from functions import simplify_fasta_name, dico

script_path = os.path.dirname(sys.argv[0])

##################
###### DEF2 ######
##################
def base_composition(seq):
      count_A=string.count(seq, "A")
      count_T=string.count(seq, "T")
      count_C=string.count(seq, "C")
      count_G=string.count(seq, "G")


      CG = count_C+count_G
      AT = count_T+count_A
      
      AG = count_A+count_G
      TC = count_T+count_C
      
      ## 1 ## Search for compositional bias in genome as marker of thermal adaptation: CG vs AT
      ratio_CG_AT=float(CG)/float(AT)
      
      ## 2 ## Search for compositional bias in genome as marker of thermal adaptation: AG vs TC
      ratio_purine_pyrimidine=float(AG)/float(TC)

      ## 3 ## Nucleotide proportion
      ln = len(seq)
      prop_A = float(count_A)/float(ln)
      prop_T = float(count_T)/float(ln)
      prop_C = float(count_C)/float(ln)
      prop_G = float(count_G)/float(ln)
      

      return(ratio_CG_AT, ratio_purine_pyrimidine, prop_A, prop_T, prop_C, prop_G)
##############################################


##################
###### DEF3 ######
##################
def aa_composition1(seq):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")
    
    

    ## 2 ## compute relative proportion
    TOTAL=count_K+count_R+count_A+count_F+count_I+count_L+count_M+count_V+count_W+count_N+count_Q+count_S+count_T+count_H+count_Y+count_C+count_D+count_E+count_P+count_G
    if (TOTAL!=0):
        ln = TOTAL
    
        prop_K=float(count_K)/float(ln)
        prop_R=float(count_R)/float(ln)
        prop_A=float(count_A)/float(ln)
        prop_F=float(count_F)/float(ln)
        prop_I=float(count_I)/float(ln)
        prop_L=float(count_L)/float(ln)
        prop_M=float(count_M)/float(ln)
        prop_V=float(count_V)/float(ln)
        prop_W=float(count_W)/float(ln)
        prop_N=float(count_N)/float(ln)
        prop_Q=float(count_Q)/float(ln)
        prop_S=float(count_S)/float(ln)
        prop_T=float(count_T)/float(ln)
        prop_H=float(count_H)/float(ln)
        prop_Y=float(count_Y)/float(ln)
        prop_C=float(count_C)/float(ln)
        prop_D=float(count_D)/float(ln)
        prop_E=float(count_E)/float(ln)
        prop_P=float(count_P)/float(ln)
        prop_G=float(count_G)/float(ln)
    else:
        prop_K=0
        prop_R=0
        prop_A=0
        prop_F=0
        prop_I=0
        prop_L=0
        prop_M=0
        prop_V=0
        prop_W=0
        prop_N=0
        prop_Q=0
        prop_S=0
        prop_T=0
        prop_H=0
        prop_Y=0
        prop_C=0
        prop_D=0
        prop_E=0
        prop_P=0
        prop_G=0



        
    return(prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G)
      
##################
###### DEF4 ######
##################
def aa_composition2(seq):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")


	
    ## 2 ## compute seq length
    TOTAL=count_K+count_R+count_A+count_F+count_I+count_L+count_M+count_V+count_W+count_N+count_Q+count_S+count_T+count_H+count_Y+count_C+count_D+count_E+count_P+count_G
    if (TOTAL!=0):

        ln = TOTAL
        ##3 Famous Hyperthermophile Prokaryotes criterias
    
        # 3.1. IVYWREL estimator  => positivelly correlated with otpimal growth
        count_IVYWREL = count_I+count_V+count_Y+count_W+count_R+count_E+count_L
        prop_IVYWREL = float(count_IVYWREL)/float(ln)
    
        # 3.2. ERK estimator (i.e. ERK vs DNQTSHA) => positivelly correlated with optimal growth temperature
        # ERK alone
        count_ERK = count_E + count_R + count_K
        prop_ERK = float(count_ERK)/float(ln)
        # DNQTSHA alone
        count_DNQTSH = count_D+count_N+count_Q+count_T+count_S+count_H
        prop_DNQTSH=float(count_DNQTSH)/float(ln)
        # ERK vs DNQTSH
        if count_DNQTSH != 0:
            ratio_ERK_vs_DNQTSH = float(count_ERK)/float(count_DNQTSH)
        else:
            ratio_ERK_vs_DNQTSH=-1
        # EK/QH estimator
        count_EK = count_E+count_K
        count_QH = count_Q+count_H

        prop_EK = float(count_EK)/float(ln)
        prop_QH = float(count_QH)/float(ln)

        if count_QH != 0:
            ratio_EK_vs_QH = float(count_EK)/float(count_QH)
        else:
            ratio_EK_vs_QH=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)

        ## 4 ## Mutationnal bias hypothesis => AT rich: favor FYMINK   // GC rich: favor GARP
        ## The mutational bias model predict a linear relationship between GARP vs FYMINK    ==> so if outliers to that, it means that the excess of GARP or FYMINK are not explained by the mutationnal bias model but by other thing ... selection!!???
        count_FYMINK=count_F+count_Y+count_M+count_I+count_N+count_K
        prop_FYMINK = float(count_FYMINK)/float(ln)

        count_GARP=count_G+count_A+count_R+count_P
        prop_GARP=float(count_GARP)/float(ln)

        ## 5 ## Hydophobicity hypothesis [should INCREASE with thermal adaptation]
        ## 5.1. AL
        count_AVLIMFYW = count_A+count_V+count_L+count_I+count_F+count_Y+count_W+count_M
        prop_AVLIMFYW=float(count_AVLIMFYW)/float(ln)
        ## 5.2. Only non-aromatic
        count_AVLIM = count_A+count_V+count_L+count_I+count_M
        prop_AVLIM=float(count_AVLIM)/float(ln)
        ## 5.3. Only aromatic (have they higher residus volume?? in such case opposite hypothesis based on residu volume, predict DECREASE for these aa in composition) 
        count_FYW = count_F+count_Y+count_W
        prop_FYW=float(count_FYW)/float(ln)

        ## 6 ## Charged hypothesis  => positivelly correlated with optimal growth temperature
        # All charged
        count_RHKDE = count_R + count_H +count_K + count_D + count_E
        prop_RHKDE = float(count_RHKDE)/float(ln)
        # Only positive
        count_RHK = count_R + count_H +count_K
        prop_RHK = float(count_RHK)/float(ln)
        # Only negative
        count_DE = count_D + count_E
        prop_DE = float(count_DE)/float(ln)
    
        ## 7 ## Neutral polar hypothesis  [should DECREASE with thermal adaptation]
        count_STNQ = count_S+count_T+count_N+count_Q
        prop_STNQ=float(count_STNQ)/float(ln)


        ## 9 ## PAYRE VS MGDS (FONTANILLAS CRITERIA)
        ## 9.1 ## Didier's criteria 1 = SMALL / BIG
        count_PAYRE = count_A+count_Y+count_P+count_R+count_E
        prop_PAYRE=float(count_PAYRE)/float(ln)
        count_MVGDS = count_V+count_M+count_S+count_G+count_D
        prop_MVGDS=float(count_MVGDS)/float(ln)
        if count_MVGDS!= 0:
            ratio_PAYRE_vs_MVGDS = float(count_PAYRE)/float(count_MVGDS)
        else:
            ratio_PAYRE_vs_MVGDS=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)

        ## 9.2 ## Didier's criteria 2 = VERY SMALL / BIG
        count_AC = count_A+count_C
        prop_AC=float(count_AC)/float(ln)

        #count_VLIM = count_V+count_L+count_I+count_M
        if count_MVGDS != 0:
            ratio_AC_vs_MVGDS = float(count_AC)/float(count_MVGDS)
        else:
            ratio_AC_vs_MVGDS=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)
    else:
        count_IVYWREL=0
        prop_IVYWREL=0
        count_ERK=0
        prop_ERK=0
        count_DNQTSH=0
        prop_DNQTSH=0
        ratio_ERK_vs_DNQTSH=0
        count_EK=0
        prop_EK=0
        count_QH=0
        prop_QH=0
        ratio_EK_vs_QH=0
        count_FYMINK=0
        prop_FYMINK=0
        count_GARP=0
        prop_GARP=0
        count_AVLIMFYW=0
        prop_AVLIMFYW=0
        count_AVLIM=0
        prop_AVLIM=0
        count_FYW=0
        prop_FYW=0
        count_STNQ=0
        prop_STNQ=0
        count_MVGDS=0
        prop_MVGDS=0
        count_PAYRE=0
        prop_PAYRE=0
        count_AC=0
        prop_AC=0
        ratio_PAYRE_vs_MVGDS=0
        ratio_AC_vs_MVGDS=0
        count_RHKDE=0
        prop_RHKDE=0
        count_RHK=0
        prop_RHK=0
        count_DE=0
        prop_DE=0

    return(count_IVYWREL,prop_IVYWREL,count_ERK,prop_ERK,count_DNQTSH,prop_DNQTSH,ratio_ERK_vs_DNQTSH,count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH,count_FYMINK,prop_FYMINK,count_GARP,prop_GARP,count_AVLIMFYW, prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW,count_STNQ, prop_STNQ, count_MVGDS,prop_MVGDS, count_PAYRE,prop_PAYRE, count_AC,prop_AC, ratio_PAYRE_vs_MVGDS, ratio_AC_vs_MVGDS, count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE)
#####################


##################
###### DEF5 ######
##################
def aa_properties(fileIN_aaProperties):
    next = fileIN_aaProperties.readline()    ## JUMP HEADERS

    bash_aa_properties={}

    while 1:
        next = fileIN_aaProperties.readline()
        if not next:
            break

        S1 = string.split(next, ",")

        aa_name = S1[1]
        S2 = string.split(aa_name, "/")
        aa_code = S2[1][:-1]

        frequencies = S1[2][:-1]
        Residue_Weight = S1[5]
        Residue_Volume = S1[6]
        Partial_specific_volume = S1[7]
        Hydration = S1[8]

        bash_aa_properties[aa_code] = [frequencies,Residue_Weight,Residue_Volume,Partial_specific_volume,Hydration]

    return(bash_aa_properties)
    
    
##################
###### DEF6 ######
##################
def sequence_properties_from_aa_properties(seq, bash_properties):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")

    TOTAL=count_K+count_R+count_A+count_F+count_I+count_L+count_M+count_V+count_W+count_N+count_Q+count_S+count_T+count_H+count_Y+count_C+count_D+count_E+count_P+count_G

    if (TOTAL!=0):


        ## 2 ## Compute properties 1: Residue Weight (Mr) (UNIT:Daltons):

        Total_Residue_Weight = count_K*float(bash_properties["K"][1]) + count_R*float(bash_properties["R"][1]) + count_A*float(bash_properties["A"][1]) + count_F*float(bash_properties["F"][1]) + count_I*float(bash_properties["I"][1]) + count_L*float(bash_properties["L"][1]) + count_M*float(bash_properties["M"][1]) + count_V*float(bash_properties["V"][1]) + count_W*float(bash_properties["W"][1]) + count_N*float(bash_properties["N"][1]) + count_Q*float(bash_properties["Q"][1]) + count_S*float(bash_properties["S"][1]) + count_T*float(bash_properties["T"][1]) + count_H*float(bash_properties["H"][1]) + count_Y*float(bash_properties["Y"][1]) + count_C*float(bash_properties["C"][1]) + count_D*float(bash_properties["D"][1]) + count_E*float(bash_properties["E"][1]) + count_P*float(bash_properties["P"][1]) + count_G*float(bash_properties["G"][1])
        Total_Residue_Volume = count_K*float(bash_properties["K"][2]) + count_R*float(bash_properties["R"][2]) + count_A*float(bash_properties["A"][2]) + count_F*float(bash_properties["F"][2]) + count_I*float(bash_properties["I"][2]) + count_L*float(bash_properties["L"][2]) + count_M*float(bash_properties["M"][2]) + count_V*float(bash_properties["V"][2]) + count_W*float(bash_properties["W"][2]) + count_N*float(bash_properties["N"][2]) + count_Q*float(bash_properties["Q"][2]) + count_S*float(bash_properties["S"][2]) + count_T*float(bash_properties["T"][2]) + count_H*float(bash_properties["H"][2]) + count_Y*float(bash_properties["Y"][2]) + count_C*float(bash_properties["C"][2]) + count_D*float(bash_properties["D"][2]) + count_E*float(bash_properties["E"][2]) + count_P*float(bash_properties["P"][2]) + count_G*float(bash_properties["G"][2])
        Total_Partial_specific_volume = count_K*float(bash_properties["K"][3]) + count_R*float(bash_properties["R"][3]) + count_A*float(bash_properties["A"][3]) + count_F*float(bash_properties["F"][3]) + count_I*float(bash_properties["I"][3]) + count_L*float(bash_properties["L"][3]) + count_M*float(bash_properties["M"][3]) + count_V*float(bash_properties["V"][3]) + count_W*float(bash_properties["W"][3]) + count_N*float(bash_properties["N"][3]) + count_Q*float(bash_properties["Q"][3]) + count_S*float(bash_properties["S"][3]) + count_T*float(bash_properties["T"][3]) + count_H*float(bash_properties["H"][3]) + count_Y*float(bash_properties["Y"][3]) + count_C*float(bash_properties["C"][3]) + count_D*float(bash_properties["D"][3]) + count_E*float(bash_properties["E"][3]) + count_P*float(bash_properties["P"][3]) + count_G*float(bash_properties["G"][3])
        Total_Hydration = count_K*float(bash_properties["K"][4]) + count_R*float(bash_properties["R"][4]) + count_A*float(bash_properties["A"][4]) + count_F*float(bash_properties["F"][4]) + count_I*float(bash_properties["I"][4]) + count_L*float(bash_properties["L"][4]) + count_M*float(bash_properties["M"][4]) + count_V*float(bash_properties["V"][4]) + count_W*float(bash_properties["W"][4]) + count_N*float(bash_properties["N"][4]) + count_Q*float(bash_properties["Q"][4]) + count_S*float(bash_properties["S"][4]) + count_T*float(bash_properties["T"][4]) + count_H*float(bash_properties["H"][4]) + count_Y*float(bash_properties["Y"][4]) + count_C*float(bash_properties["C"][4]) + count_D*float(bash_properties["D"][4]) + count_E*float(bash_properties["E"][4]) + count_P*float(bash_properties["P"][4]) + count_G*float(bash_properties["G"][4])
    else:
        Total_Residue_Weight=0
        Total_Residue_Volume=0
        Total_Partial_specific_volume=0
        Total_Hydration=0
    
    return(Total_Residue_Weight,Total_Residue_Volume,Total_Partial_specific_volume,Total_Hydration)

########################################################



###################
### RUN RUN RUN ###
###################


##Create specific folders
Path_IN_loci_NUC = "./IN_AA"
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


## 2 ## PathIN
fileIN_properties = open("amino_acid_properties.csv", "r")
Path_IN_loci_AA = "./IN_AA"
#Path_IN_loci_AA = "02_CDS_No_Missing_Data_aa_CDS_withM"
Lloci_AA = os.listdir(Path_IN_loci_AA)

## 3 ## PathOUT

## 3.1 ## PROT composition
fileOUT_PROT_ALL=open("./OUT/prot_compositions_All_AA.csv","w")
fileOUT_PROT_ALL.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_PROT_ALL.write("%s_prop_K,%s_prop_R,%s_prop_A,%s_prop_F,%s_prop_I,%s_prop_L,%s_prop_M,%s_prop_V,%s_prop_W,%s_prop_N,%s_prop_Q,%s_prop_S,%s_prop_T,%s_prop_H,%s_prop_Y,%s_prop_C,%s_prop_D,%s_prop_E,%s_prop_P,%s_prop_G," %(taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_PROT_ALL.write("%s_prop_K,%s_prop_R,%s_prop_A,%s_prop_F,%s_prop_I,%s_prop_L,%s_prop_M,%s_prop_V,%s_prop_W,%s_prop_N,%s_prop_Q,%s_prop_S,%s_prop_T,%s_prop_H,%s_prop_Y,%s_prop_C,%s_prop_D,%s_prop_E,%s_prop_P,%s_prop_G" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_PROT_ALL.write("\n")

## 3.2 ## PROT IVYWREL
fileOUT_IVYWREL=open("./OUT/IVYWREL.csv","w")
fileOUT_IVYWREL.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_IVYWREL.write("%s_count_IVYWREL,%s_prop_IVYWREL," %(taxa,taxa))
fileOUT_IVYWREL.write("%s_count_IVYWREL,%s_prop_IVYWREL" %(LT[-1],LT[-1]))
fileOUT_IVYWREL.write("\n")

## 3.3 ## PROT ERK_DNQTSHA
fileOUT_ERK_DNQTSH=open("./OUT/ERK_DNQTSH.csv","w")
fileOUT_ERK_DNQTSH.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_ERK_DNQTSH.write("%s_count_ERK,%s_prop_ERK,%s_count_DNQTSH,%s_prop_DNQTSH,%s_ratio_ERK_vs_DNQTSH," %(taxa,taxa,taxa,taxa,taxa))
fileOUT_ERK_DNQTSH.write("%s_count_ERK,%s_prop_ERK,%s_count_DNQTSH,%s_prop_DNQTSH,%s_ratio_ERK_vs_DNQTSH" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_ERK_DNQTSH.write("\n")

## 3.4 ## PROT EK_QH
fileOUT_EK_QH=open("./OUT/EK_QH.csv","w")
fileOUT_EK_QH.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_EK_QH.write("%s_count_EK,%s_prop_EK,%s_count_QH,%s_prop_QH,%s_ratio_EK_vs_QH," %(taxa,taxa,taxa,taxa,taxa))
fileOUT_EK_QH.write("%s_count_EK,%s_prop_EK,%s_count_QH,%s_prop_QH,%s_ratio_EK_vs_QH" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_EK_QH.write("\n")

## 3.5 ## PROT FYMINK_GARP
fileOUT_FYMINK_GARP=open("./OUT/FYMINK_GARP.csv","w")
fileOUT_FYMINK_GARP.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_FYMINK_GARP.write("%s_count_FYMINK,%s_prop_FYMINK,%s_count_GARP,%s_prop_GARP," %(taxa,taxa,taxa,taxa))
fileOUT_FYMINK_GARP.write("%s_count_FYMINK,%s_prop_FYMINK,%s_count_GARP,%s_prop_GARP" %(LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_FYMINK_GARP.write("\n")

## 3.6 ## PROT AVLIMFYW
fileOUT_AVLIMFYW=open("./OUT/AVLIMFYW.csv","w")
fileOUT_AVLIMFYW.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_AVLIMFYW.write("%s_count_AVLIMFYW,%s_prop_AVLIMFYW,%s_count_AVLIM,%s_prop_AVLIM,%s_count_FYW,%s_prop_FYW," %(taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_AVLIMFYW.write("%s_count_AVLIMFYW,%s_prop_AVLIMFYW,%s_count_AVLIM,%s_prop_AVLIM,%s_count_FYW,%s_prop_FYW" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_AVLIMFYW.write("\n")

## 3.7 ## PROT STNQ
fileOUT_STNQ=open("./OUT/STNQ.csv","w")
fileOUT_STNQ.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_STNQ.write("%s_count_STNQ,%s_prop_STNQ," %(taxa,taxa))
fileOUT_STNQ.write("%s_count_STNQ,%s_prop_STNQ" %(LT[-1],LT[-1]))
fileOUT_STNQ.write("\n")

## 3.8 ## PROT RHKDE
fileOUT_RHKDE=open("./OUT/RHKDE.csv","w")
fileOUT_RHKDE.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_RHKDE.write("%s_count_RHKDE,%s_prop_RHKDE,%s_count_RHK,%s_prop_RHK,%s_count_DE,%s_prop_DE," %(taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_RHKDE.write("%s_count_RHKDE,%s_prop_RHKDE,%s_count_RHK,%s_prop_RHK,%s_count_DE,%s_prop_DE" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_RHKDE.write("\n")

## 3.9 ## PROT DIDER CRITERIA
fileOUT_PAYRE=open("./OUT/PAYRE-MVGDS.csv","w")
fileOUT_PAYRE.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_PAYRE.write("%s_count_PAYRE,%s_prop_PAYRE,%s_count_AC,%s_prop_AC,%s_count_MVGDS,%s_prop_MVGDS,%s_ratio_PAYRE_vs_MVGDS,%s_ratio_AC_vs_MVGDS," %(taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_PAYRE.write("%s_count_PAYRE,%s_prop_PAYRE,%s_count_AC,%s_prop_AC,%s_count_MVGDS,%s_prop_MVGDS,%s_ratio_PAYRE_vs_MVGDS,%s_ratio_AC_vs_MVGDS" %(LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1],LT[-1]))
fileOUT_PAYRE.write("\n")

## 3.10 ## PROT Total residue weight
fileOUT_TotalResidueWeight=open("./OUT/TotalResidueWeight.csv","w")
fileOUT_TotalResidueWeight.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_TotalResidueWeight.write("%s_Total_Residue_Weight," %taxa)
fileOUT_TotalResidueWeight.write("%s_Total_Residue_Weight" %LT[-1])
fileOUT_TotalResidueWeight.write("\n")

## 3.11 ## PROT Total residue volume
fileOUT_TotalResidueVolume=open("./OUT/TotalResidueVolume.csv","w")
fileOUT_TotalResidueVolume.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_TotalResidueVolume.write("%s_Total_Residue_Volume," %taxa)
fileOUT_TotalResidueVolume.write("%s_Total_Residue_Volume" %LT[-1])
fileOUT_TotalResidueVolume.write("\n")

## 3.12 ## PROT Total partial specific volume
fileOUT_TotalPartialSpecificVolume=open("./OUT/TotalPartialSpecificVolume.csv","w")
fileOUT_TotalPartialSpecificVolume.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_TotalPartialSpecificVolume.write("%s_Total_Partial_Specific_Volume," %taxa)
fileOUT_TotalPartialSpecificVolume.write("%s_Total_Partial_Specific_Volume" %LT[-1])
fileOUT_TotalPartialSpecificVolume.write("\n")

## 3.13 ## PROT Total hydratation
fileOUT_TotalHydratation=open("./OUT/TotalHydratation.csv","w")
fileOUT_TotalHydratation.write("LOCUS,")
for taxa in LT[0:-1]:
    fileOUT_TotalHydratation.write("%s_Total_Hydratation," %taxa)
fileOUT_TotalHydratation.write("%s_Total_Hydratation" %LT[-1])
fileOUT_TotalHydratation.write("\n")

#####################
## 4 ## Process Loci
#####################
bash_aa_properties = aa_properties(fileIN_properties)

for locus in Lloci_AA:
    if not os.path.isfile("%s/%s" %(Path_IN_loci_AA, locus)) :
        pass
    else :
        print locus
        path_locus = "%s/%s" %(Path_IN_loci_AA, locus)
        bash = dico(path_locus,LT)

        #print bash
        
        fileOUT_PROT_ALL.write("%s," %locus)
        fileOUT_IVYWREL.write("%s," %locus)
        fileOUT_ERK_DNQTSH.write("%s," %locus)
        fileOUT_EK_QH.write("%s," %locus)
        fileOUT_FYMINK_GARP.write("%s," %locus)
        fileOUT_AVLIMFYW.write("%s," %locus)
        fileOUT_STNQ.write("%s," %locus)
        fileOUT_RHKDE.write("%s," %locus)
        fileOUT_PAYRE.write("%s," %locus)
        fileOUT_TotalResidueWeight.write("%s," %locus)
        fileOUT_TotalResidueVolume.write("%s," %locus)
        fileOUT_TotalPartialSpecificVolume.write("%s," %locus)
        fileOUT_TotalHydratation.write("%s," %locus)
        
        for taxa in LT[0:-1]:
            if taxa in bash.keys():
                seq = bash[taxa]
                prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G = aa_composition1(seq)   ### DEF3 ###
                count_IVYWREL,prop_IVYWREL,count_ERK,prop_ERK,count_DNQTSH,prop_DNQTSH,ratio_ERK_vs_DNQTSH,count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH,count_FYMINK,prop_FYMINK,count_GARP,prop_GARP,count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW,count_STNQ,prop_STNQ, count_MVGDS,prop_MVGDS, count_PAYRE,prop_PAYRE, count_AC,prop_AC, ratio_PAYRE_vs_MVGDS, ratio_AC_vs_MVGDS,count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE = aa_composition2(seq)   ### DEF4 ###
                Total_Residue_Weight,Total_Residue_Volume,Total_Partial_Specific_Volume,Total_Hydration = sequence_properties_from_aa_properties(seq, bash_aa_properties)   ### DEF6 ###
            
                fileOUT_PROT_ALL.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G))
                fileOUT_IVYWREL.write("%.5f,%.5f," %(count_IVYWREL, prop_IVYWREL))
                fileOUT_ERK_DNQTSH.write("%.5f,%.5f,%.5f,%.5f,%.5f," %(count_ERK,prop_ERK,count_DNQTSH,prop_DNQTSH,ratio_ERK_vs_DNQTSH))
                fileOUT_EK_QH.write("%.5f,%.5f,%.5f,%.5f,%.5f," %(count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH))
                fileOUT_FYMINK_GARP.write("%.5f,%.5f,%.5f,%.5f," %(count_FYMINK,prop_FYMINK,count_GARP,prop_GARP))
                fileOUT_AVLIMFYW.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW))
                fileOUT_STNQ.write("%.5f,%.5f," %(count_STNQ,prop_STNQ))
                fileOUT_RHKDE.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,"%(count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE))
                fileOUT_PAYRE.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(count_PAYRE,prop_PAYRE,count_AC,prop_AC,count_MVGDS,prop_MVGDS,ratio_PAYRE_vs_MVGDS,ratio_AC_vs_MVGDS))
                fileOUT_TotalResidueWeight.write("%.5f," %Total_Residue_Weight)
                fileOUT_TotalResidueVolume.write("%.5f," %Total_Residue_Volume)
                fileOUT_TotalPartialSpecificVolume.write("%.5f," %(Total_Partial_Specific_Volume))
                fileOUT_TotalHydratation.write("%.5f," % Total_Hydration)
            else:
                fileOUT_PROT_ALL.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"))
                fileOUT_IVYWREL.write("%s,%s," %("NA", "NA"))
                fileOUT_ERK_DNQTSH.write("%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA"))
                fileOUT_EK_QH.write("%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA"))
                fileOUT_FYMINK_GARP.write("%s,%s,%s,%s," %("NA","NA","NA","NA"))
                fileOUT_AVLIMFYW.write("%s,%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA","NA"))
                fileOUT_STNQ.write("%s,%s," %("NA","NA"))
                fileOUT_RHKDE.write("%s,%s,%s,%s,%s,%s,"%("NA","NA","NA","NA","NA","NA"))
                fileOUT_PAYRE.write("%s,%s,%s,%s,%s,%s,%s,%s," %("NA","NA","NA","NA","NA","NA","NA","NA"))
                fileOUT_TotalResidueWeight.write("%s," %"NA")
                fileOUT_TotalResidueVolume.write("%s," %"NA")
                fileOUT_TotalPartialSpecificVolume.write("%s," %"NA")
                fileOUT_TotalHydratation.write("%s," %"NA")

        if LT[-1] in bash.keys():
            seq = bash[LT[-1]]            
            prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G = aa_composition1(seq)   ### DEF3 ###
            count_IVYWREL,prop_IVYWREL,count_ERK,prop_ERK,count_DNQTSH,prop_DNQTSH,ratio_ERK_vs_DNQTSH,count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH,count_FYMINK,prop_FYMINK,count_GARP,prop_GARP,count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW,count_STNQ,prop_STNQ, count_MVGDS,prop_MVGDS, count_PAYRE,prop_PAYRE, count_AC,prop_AC, ratio_PAYRE_vs_MVGDS, ratio_AC_vs_MVGDS,count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE = aa_composition2(seq)   ### DEF4 ###
            Total_Residue_Weight,Total_Residue_Volume,Total_Partial_Specific_Volume,Total_Hydration = sequence_properties_from_aa_properties(seq, bash_aa_properties)   ### DEF6 ###
        
            fileOUT_PROT_ALL.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f" %(prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G))
            fileOUT_IVYWREL.write("%.5f,%.5f" %(count_IVYWREL, prop_IVYWREL))
            fileOUT_ERK_DNQTSH.write("%.5f,%.5f,%.5f,%.5f,%.5f" %(count_ERK,prop_ERK,count_DNQTSH,prop_DNQTSH,ratio_ERK_vs_DNQTSH))
            fileOUT_EK_QH.write("%.5f,%.5f,%.5f,%.5f,%.5f" %(count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH))
            fileOUT_FYMINK_GARP.write("%.5f,%.5f,%.5f,%.5f" %(count_FYMINK,prop_FYMINK,count_GARP,prop_GARP))
            fileOUT_AVLIMFYW.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f" %(count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW))
            fileOUT_STNQ.write("%.5f,%.5f" %(count_STNQ,prop_STNQ))
            fileOUT_RHKDE.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f"%(count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE))
            fileOUT_PAYRE.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f" %(count_PAYRE,prop_PAYRE,count_AC,prop_AC,count_MVGDS,prop_MVGDS,ratio_PAYRE_vs_MVGDS,ratio_AC_vs_MVGDS))
            fileOUT_TotalResidueWeight.write("%.5f" %Total_Residue_Weight)
            fileOUT_TotalResidueVolume.write("%.5f" %Total_Residue_Volume)
            fileOUT_TotalPartialSpecificVolume.write("%.5f" %(Total_Partial_Specific_Volume))
            fileOUT_TotalHydratation.write("%.5f" % Total_Hydration)
        else:
            fileOUT_PROT_ALL.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"))
            fileOUT_IVYWREL.write("%s,%s" %("NA", "NA"))
            fileOUT_ERK_DNQTSH.write("%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA"))
            fileOUT_EK_QH.write("%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA"))
            fileOUT_FYMINK_GARP.write("%s,%s,%s,%s" %("NA","NA","NA","NA"))
            fileOUT_AVLIMFYW.write("%s,%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA","NA"))
            fileOUT_STNQ.write("%s,%s" %("NA","NA"))
            fileOUT_RHKDE.write("%s,%s,%s,%s,%s,%s"%("NA","NA","NA","NA","NA","NA"))
            fileOUT_PAYRE.write("%s,%s,%s,%s,%s,%s,%s,%s" %("NA","NA","NA","NA","NA","NA","NA","NA"))
            fileOUT_TotalResidueWeight.write("%s" %"NA")
            fileOUT_TotalResidueVolume.write("%s" %"NA")
            fileOUT_TotalPartialSpecificVolume.write("%s" %"NA")
            fileOUT_TotalHydratation.write("%s" %"NA")
    
    ## END LINE
    fileOUT_PROT_ALL.write("\n")
    fileOUT_IVYWREL.write("\n")
    fileOUT_ERK_DNQTSH.write("\n")
    fileOUT_EK_QH.write("\n")
    fileOUT_FYMINK_GARP.write("\n")
    fileOUT_AVLIMFYW.write("\n")
    fileOUT_STNQ.write("\n")
    fileOUT_RHKDE.write("\n")    
    fileOUT_PAYRE.write("\n")
    fileOUT_TotalResidueWeight.write("\n")
    fileOUT_TotalResidueVolume.write("\n")
    fileOUT_TotalPartialSpecificVolume.write("\n")
    fileOUT_TotalHydratation.write("\n")

