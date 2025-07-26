# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 19:24:26 2022

@author: selcuk.1
"""


def dict_to_fasta(MSA_dict,output):
    file=open(output,"w")
    for header in MSA_dict:
        file.write(header+"\n")
        file.write(MSA_dict[header]+"\n")
    file.close()
    
def get_cluster_dict(cluster_file):
    cluster_id_d={}
    cluster_read=open(cluster_file,"r")
    cl_count=-1
    for line in cluster_read:
        if line[0]==">":
            cl_count+=1
            key="Cl"+str(cl_count)
            cluster_id_d[key]=[]
        else:
            line_list=line.split("|")
            ID=line_list[1]
            cluster_id_d[key].append(ID)
    cluster_read.close()
    return cluster_id_d            
    

def sort_list(list1, list2):
 
    zipped_pairs = zip(list2, list1)
 
    z = [x for _, x in sorted(zipped_pairs)]
     
    return z
    
def get_cluster_sizes(cluster_dict):
    length_list=[]
    cluster_name_list=cluster_dict.keys()
    for key in cluster_dict:
        length_cluster=len(cluster_dict[key])
        length_list.append(length_cluster)
    sorted_keys=sort_list(cluster_name_list,length_list)
    return sorted_keys
        
import os
def get_reps(receptor_list,topn):
    output_dict={}
    for receptor in receptor_list:
        for i in range(0,8):
            class_list=["_classA","_classB","_classC","_classF","_classT","_classU","_outgroup","_classDandE"]
            class_id=class_list[i]
            date="26-10-2022"
            if class_id=="_classDandE":
                directory=r"/cta/users/bselcuk/GPCRA/GPCR_atlas_26-10-2022/{}_{}{}/{}_{}_blast_reps_0.65.fasta".format(receptor,date,class_id,receptor,date)
            else:
                directory=r"/cta/users/bselcuk/GPCRA/GPCR_atlas_26-10-2022/{}_{}{}/{}_{}_subtree_reps_0.65.fasta".format(receptor,date,class_id,receptor,date)
            file_exists= os.path.isfile(directory)
            if file_exists:
                if class_id=="_classDandE":
                    cluster_dict=get_cluster_dict(r"/cta/users/bselcuk/GPCRA/GPCR_atlas_26-10-2022/{}_{}{}/{}_{}_blast_reps_0.65.fasta.clstr".format(receptor,date,class_id,receptor,date))
                else:
                    cluster_dict=get_cluster_dict(r"/cta/users/bselcuk/GPCRA/GPCR_atlas_26-10-2022/{}_{}{}/{}_{}_subtree_reps_0.65.fasta.clstr".format(receptor,date,class_id,receptor,date))
                break
        sorted_clusters=get_cluster_sizes(cluster_dict)
        rep_fasta=open(directory,"r")
        seq_count=0
        cluster_names=sorted_clusters[(topn*-1):]
        for line in rep_fasta:
            line=line.strip()
            line_list=line.split("|")
            if line[0]==">":
                write_seq=0
                ID=line_list[1]
                for name in cluster_names:
                    IDs=cluster_dict[name]
                    if ID in IDs:
                        write_seq=1
                        key=line+"_{}".format(receptor)
                        output_dict[key]=""
                        break
            elif write_seq==1:
                output_dict[key]+=line
    return output_dict

def human_seq_get(class_id,receptor_list,sequence_dict):
    date="26-10-2022"
    for receptor in receptor_list:
        sequence=""
        directory=directory=r"/cta/users/bselcuk/GPCRA/GPCR_atlas_26-10-2022/{}_{}{}/{}_{}_subtree.fasta".format(receptor,date,class_id,receptor,date)
        fasta_file=open(directory,"r")
        write=0
        for line in fasta_file:
            if "9606" in line and receptor in line:
                header=line.strip("\n")
                write=1
            elif write==1 and "|" in line:
                fasta_file.close()
                break
            elif write==1:
                sequence+=line.strip("\n")
        sequence_dict[header]=sequence
    return sequence_dict
    
all_human_classA=['5HT1A', 'OXGR1', 'CNR1', 'NK2R', 'CCR1', '5HT1E', 'GPR4', 'OPSR', 'GPR62', 'PKR2', 'MRGRE', 'MTLR', 'OPSG2', 'GRPR', 'XCR1', 'SSR5', 'GPR3', 'CCR7', 'NPY4R', 'MCHR2', 'P2Y10', 'GPR83', 'GPR6', 'RXFP1', 'NPSR1', 'G37L1', '5HT2B', 'TSHR', 'DRD2', '5HT7R', 'HRH4', 'RL3R1', 'NMUR1', 'HCAR2', 'C3AR', 'OPSB', 'OX2R', 'GPR45', 'P2RY6', 'CML1', 'GPR27', 'OXER1', 'GPR87', 'DRD1', 'TAAR5', 'TAAR9', 'PE2R2', 'HRH2', 'NPBW2', 'ADRB2', 'OPRX', 'PI2R', 'GP139', 'MC5R', 'PAR1', 'GPR84', 'CCR2', 'V1AR', 'NK3R', 'GP119', 'TAAR1', 'GP146', 'ADA1D', 'FFAR3', 'GPR88', 'GPR82', 'OPRD', 'OPSG', '5HT2C', 'CXCR1', 'CCR6', 'P2Y14', 'ACTHR', 'CCKAR', 'MSHR', 'MTR1B', 'FFAR1', 'GP142', 'PAR3', 'FFAR2', 'V1BR', 'OPSD', '5HT4R', 'FFAR4', 'PSYR', 'UR2R', 'PTAFR', 'GPR1', 'GPR63', 'ADA2A', 'HRH3', 'MRGRD', 'ADA2B', 'NMBR', 'CX3C1', 'P2Y13', 'GP151', 'C5AR1', 'NTR2', 'GPR12', 'PD2R2', 'TAAR3', 'GP101', 'GASR', 'LPAR5', 'AGTR2', 'GHSR', 'GPR26', 'OXYR', 'LPAR6', 'TAAR6', 'ACKR4', 'LGR6', 'OPRK', 'AA2BR', 'CCR10', 'ACM4', 'SSR1', 'NMUR2', 'GPR78', 'GP141', 'GNRHR', 'LGR4', 'SSR2', 'OPN3', 'MC3R', 'NPY2R', 'LPAR4', 'AA1R', 'GPR21', 'AA3R', 'HCAR1', 'NPY6R', 'GPR33', 'NPY1R', 'S1PR4', 'MCHR1', 'EDNRB', 'NPFF2', 'RL3R2', 'DRD5', 'ADRB1', 'FPR2', 'CXCR3', 'CCR8', 'CCRL2', 'NK1R', 'FPR1', 'GPR15', 'G32P1', 'ACKR2', 'TAAR8', 'ACKR3', 'NPY42', 'P2RY8', 'MRGX3', 'QRFPR', 'GPR61', 'GPR19', 'GPR22', '5HT6R', 'CXCR4', 'EDNRA', 'MTR1L', 'SSR3', 'GALR1', 'GPR55', 'MRGX2', 'ACM5', 'GPR39', 'MRGX4', 'GPR35', 'S1PR3', 'OPN4', '5HT5A', 'MRGX1', 'CCR4', 'TAAR2', 'TRFR', 'CLTR2', 'ACM1', 'GP135', 'GPR52', 'P2Y12', 'GPR85', 'PAR4', 'HRH1', 'LPAR1', 'C5AR2', 'S1PR5', 'NPY5R', 'NPBW1', 'ACM3', 'GPER1', 'CCR9', 'ADA1B', 'P2RY4', 'RGR', 'OX1R', 'PF2R', 'ADA1A', 'GPR34', 'PE2R3', 'GP183', 'GP176', 'CXCR6', 'RXFP2', 'GP162', 'ADRB3', 'PE2R1', 'GPR18', 'KISSR', 'GP174', 'FSHR', 'TA2R', 'PRLHR', 'GP182', 'BRS3', 'PKR1', 'GP161', 'GPR20', 'GALR3', 'GP150', 'CNR2', 'P2RY2', 'CCR5', 'MRGRG', 'V2R', 'GP153', '5HT1D', 'MAS', 'GPR32', 'GNRR2', 'MRGRF', 'BKRB2', 'LT4R2', 'HCAR3', 'GP132', 'LT4R1', '5HT2A', 'PE2R4', 'GP171', 'LPAR3', 'GPR17', 'S1PR1', 'GALR2', 'PAR2', 'GPBAR', 'LPAR2', 'FPR3', 'GPR75', 'P2Y11', 'OPSX', 'ACM2', 'AGTR1', 'APJ', 'CXCR2', 'GP160', 'OGR1', 'GPR42', 'AA2AR', 'PD2R', 'SSR4', 'NTR1', 'S1PR2', 'CXCR5', 'GPR37', 'SUCR1', 'MC4R', 'P2RY1', 'GPR31', 'DRD4', 'GPR25', 'GP173', 'MAS1L', 'OPRM', 'OPN5', 'OPSG3', 'ADA2C', 'LSHR', 'NPFF1', 'LGR5', '5HT1F', 'CCR3', 'MTR1A', 'ACKR1', 'BKRB1', '5HT1B', 'GP149', 'CLTR1', 'DRD3', 'GP152',"GP148"]

all_classC="CASR GABR2 GRM3 GRM2 GRM6 GRM7 GRM4 GRM8 GABR1 GRM5 GRM1 GPC6A GPC5B GP179 GP158 TS1R2 TS1R1 GPC5C GP156 TS1R3 GPC5D RAI3".split(" ")

all_classB="PACR VIPR1 SCTR PTH1R PTH2R VIPR2 CALCR GIPR GHRHR GLR GLP2R CALRL GLP1R AGRL3 CRFR2 AGRE5 AGRL2 CRFR1 AGRL1 AGRE1 AGRB1 AGRG2 AGRB3 AGRE3 AGRL4 CELR3 AGRD1 CELR1 AGRE4 AGRE2 AGRG4 AGRG6 CELR2 AGRD2 AGRF1 AGRF5 AGRG3 AGRF3 AGRF4 AGRG5 AGRG1 AGRF2 AGRA1 AGRA3 AGRV1 AGRA2 AGRB2 AGRG7".split(" ")
all_classF="FZD1 FZD2 FZD3 FZD4 FZD5 FZD6 FZD7 FZD8 FZD9 FZD10 SMO".split(" ")
all_classT="T2R43 T2R40 T2R41 T2R42 TA2R1 TA2R3 TA2R4 TA2R5 TA2R7 TA2R8 TA2R9 T2R10 T2R13 T2R14 T2R16 T2R19 T2R20 T2R30 T2R31 T2R38 T2R39 T2R45 T2R46 T2R50 T2R60".split(" ")
olfactories=["O52B2","O52H1","O52B4","O52B6","O52Z1","O52D1","O52E1","O52E2","O52E5","O52E4","O52E8","O52E6","O52J3","O52P1","O52R1","O52L1","O52L2","O52W1","O52M1","O52K1","O52K2","O51F2","O51F1","A0A3B3IT45","O51T1","O51I1","O51I2","O51J1","O51M1","O51Q1","O51L1","O51H1","O51B2","O51B6","O51B4","O51B5","O51V1","O51G2","O51G1","O51S1","O51A2","O51A4","O51A7","O51E2","O51D1","O51E1","O52A5","O52A1","O52A4","O52N4","O52N1","O52N2","O52N5","O52I2","O52I1","O56A3","O56A4","O56A5","O56A1","O56B1","O56B2","O56B4","O2AT4","O11H2","O11H7","O11H4","O11H6","O11G2","O11A1","OR6N1","OR6N2","OR6K2","OR6K3","OR6K6","OR6F1","OR9A2","OR9A1","OR9A4","OR6V1","OR6S1","OR6M1","OR6J1","OR6T1","O6C75","OR6C3","OR6C1","O6C76","O6C65","O6C74","OR6C2","O6C68","O6C70","OR6C6","OR6C4","O2AP1","OR6X1","OR6B2","OR6B3","OR6A2","OR6P1","OR6Y1","OR6B1","OR6Q1","O11L1","O2T10","OR2T4","O2T29","OR2T5","OR2T3","O2T27","OR2T7","OR2T1","OR2T6","O2T35","OR2T2","O2T11","O2T33","OR2T8","O2T33","OR2M3","OR2M7","OR2M5","OR2M2","OR2M4","OR2V1","OR2V2","O2AK2","OR2L2","OR2L5","OR2L8","OR2L3","OR2LD","O2AJ1","O2AE1","OR2Z1","O2AG1","O2AG2","OR2B2","OR2B6","OR2B3","OR2BB","OR2Y1","OR2I1","OR2B8","OR2J3","OR2J1","OR2J2","OR2W1","OR2W3","OR2W5","OR2W6","OR2C1","OR2H2","OR2H1","OR2G3","OR2C3","OR2G2","OR2G6","A0A126GWI2","OR2A4","OR2A7","O2A25","O2A14","A0A4W9AIG4","O2A12","A0A126GWB0","OR2A2","OR2A5","A0A2R8YEH3","A0A2R8YEH3","OR2D3","OR2D2","OR2F1","OR2F2","O13H1","O13F1","OR2K2","O13D1","O13C9","O13C5","O13C2","O13C6","O13C7","OR2S1","O13C8","O13C3","O13C4","O13J1","O10AD","O10A3","O10A6","O10A5","O10A2","O10A4","O10A7","O10C1","O10P1","O10AG","O10H1","O10H5","O10H2","O10H4","O10H3","O10X1","O10J4","O10J3","O10J6","O10J5","O10J1","O10K2","O10K1","O10T2","O10R2","O10Z1","O10W1","O10Q1","O10V1","O10AC","O14CZ","O14I1","O14K1","O14A2","O14AG","O14J1","O14L1","A0A2R8YED5","OR8S1","OR5V1","OR3A2","OR3A1","OR3A3","OR3A4","O13G1","O13A1","OR4E2","OR4E1","OR4DA","OR4DB","OR4D9","OR4D6","OR4D1","OR4D2","OR4D5","OR4Q3","OR4L1","OR4KD","OR4K3","OR4K2","OR4KE","OR4KH","OR4K5","OR4KF","OR4K1","OR4F6","O4F15","OR4F3","O4F21","OR4F4","OR4F4","OR4Q2","OR4N5","A0A0G2JNH3","A0A0G2JNH3","OR4N2","ORM2B","OR4M1","OR4M2","OR4CB","OR4CG","OR4C3","OR4CF","OR4A8","OR4A5","O4A16","O4A47","OR4A4","O4A15","OR4C6","OR4C5","O4C45","O4C46","OR4CD","OR4CC","OR4P4","OR4S2","OR4S1","OR4X2","OR4X1","OR4B1","O10G4","O10G7","O10G8","O10G9","O10G6","O10G3","O10G2","O10S1","O10D4","O10D3","O12D3","O12D2","O12D1","OR5P2","OR5P3","OR5G3","O5AU1","OR9G4","OR9G1","OR9G9","OR9Q2","OR9Q1","OR9I1","O5AK2","O5AK3","OR5B2","OR5B3","OR5BH","OR5BC","OR5BL","OR8D4","OR8A1","OR8D2","OR8D1","OR83P","O8G2P","OR8G5","OR8G1","OR8B4","OR8BC","OR8B8","OR8B2","OR8B3","OR5F1","OR8H2","OR8H3","OR8H1","OR8I2","OR8K5","OR8K3","OR8K1","OR8J3","OR8J1","OR8J2","OR8U1","OR8U8","OR8U9","OR5R1","O5AL1","OR5M9","OR5M3","OR5M8","OR5M1","OR5MA","OR5MB","O5AP2","OR5J2","OR9K2","OR5L2","OR5L1","OR5DI","OR5DG","OR5DE","OR5DD","A0A2R8Y4L6","O5AR1","OR5C1","OR5W2","OR5I1","O5AS1","O5H14","OR5H6","OR5H1","O5H15","OR5H8","OR5H2","O5AC1","O5AC2","OR5K1","OR5K2","OR5K3","OR5K4","OR5T2","OR5T1","OR5T3","O5AN1","OR5A2","OR5A1","A0A286YEU6","OR1Q1","OR1B1","OR1K1","OR1L6","OR1L4","OR1L1","OR1L3","OR1L8","OR1A2","OR1A1","OR1F1","OR1F2","OR1E3","OR1E2","OR1E1","OR1J1","OR1J2","OR1J4","OR1S2","OR1S1","OR1M1","OR1P1","OR1D5","OR1D4","OR1D2","OR1N2","OR1N1","OR1C1","OR1FC","OR1I1","OR1G1","OR7G1","OR7G2","OR7G3","OR7D2","OR7D4","OR7A5","OR7A2","OR7AA","OR7AH","OR7C1","OR7C2","O7E24","OR4N4","ORN4C","O4F21","OR4F3","OR4F5","O4F17","O13C7","A0A2R8YEG4","A0A2R8YEV3","OR2A1","OR2A2","A0A126GWB0","O2A25","OR2A7","OR2A4","A0A126GWI2","A0A4W9AIG4","O2A14","OR2T5","O2T29","O2T34","OR2T7","O2T27","OR2T2","O2T35","O2T12","O11H1","O11HC"]

receptor_list=all_human_classA
receptor_list=all_classC
receptor_list=all_classB
receptor_list=all_classF
receptor_list=all_classT
receptor_list=olfactories
class_id="_classA"
class_id="_classC"
class_id="_classB"
class_id="_classF"
class_id="_classT"
class_id="_classA"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classA_top5_reps_whuman_corrected.fasta"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classC_top5_reps_whuman.fasta"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classB_top5_reps_whuman.fasta"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classF_top5_reps_whuman.fasta"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classT_top5_reps_whuman.fasta"
output_dir=r"/cta/users/bselcuk/GPCRA/class_rep_alignments/classOlf_top5_reps_whuman.fasta"
reps_d=get_reps(receptor_list,5)
humans_and_reps=human_seq_get(class_id,receptor_list,reps_d)
dict_to_fasta(humans_and_reps,output_dir)





