# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:48:06 2022

@author: selcuk.1
"""
import re
def iqtree_best_model_get(log_file):
    log_read=open(log_file,"r",encoding="utf8")
    for line in log_read:
        #print(line)
        text="Bayesian Information Criterion:         "
        if line[0:len(text)]==text:
            #print("YAYY!")
            model=line[len(text):-1].strip()
    log_read.close()
    return model


import argparse
ap = argparse.ArgumentParser(description="Script for obtaining the best substitution model from IQ-Tree2 Modelfinder log file.")
ap.add_argument("--logfile", required=True,help="IQ-Tree2 modelfinder log file")
ap.add_argument("--modelout",required=True,help="Output text file including best fitting model")
args = vars(ap.parse_args())

log=args["logfile"]
out_text=args["modelout"]

final_output=open(out_text,"w")
final_output.write(iqtree_best_model_get(log))
final_output.close()
