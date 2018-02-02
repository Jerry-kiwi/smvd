# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from Bio import Entrez,SeqUtils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_nucleotide
from Bio.SeqUtils import MeltingTemp
import time
import itertools
from itertools import count # izip for maximum efficiency
import multiprocessing
from multiprocessing import Pool as ThreadPool
import os
import subprocess
import shutil

from sqlalchemy import create_engine
import sqlite3
import math

xlsx = pd.ExcelFile('SRC_DATA.xlsx')
df1 = xlsx.parse('HUMAN')
df1a = xlsx.parse('VIRUS')

blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_04b.db') 

GENE_LIST=np.unique(df1[2].values)


OUTPUT1=[]
OUTPUT2=[]
for G in GENE_LIST:
    block=df1[df1[2]==G]

    for index,row in block.iterrows():
        M=row[5]

        
        df2 = pd.read_sql_query('SELECT DISTINCT("SEQ"),"Tm","Gc","GENE","MRNA","TYPE","RATE","LEN","POS" from "main"."table" WHERE "GENE"=' + str(G) + ' AND "MRNA"="' + str(M) + '"  AND "RATE" >60 ORDER BY "RATE" DESC,"LEN" DESC', blast_db)
        
        if df2.shape[0]>0:
            df2=df2.sort_values(by=['POS','LEN'],ascending=[False,False])
            ZA=df2.iloc[0]
            S_LOW=ZA['SEQ']
            R1=ZA['RATE']       
            OUTPUT1.append([row[0],G,row[3],row[1],M,row[3]+'-'+M+'-Z1',S_LOW,len(S_LOW),R1])
            if df2.shape[0]>1: 
                df2=df2.sort_values(by=['POS','LEN'],ascending=[True,False])
                ZA=df2.iloc[0]
                S_LOW=ZA['SEQ']
                R1=ZA['RATE']       
                OUTPUT1.append([row[0],G,row[3],row[1],M,row[3]+'-'+M+'-Z2',S_LOW,len(S_LOW),R1])
            else:
                OUTPUT1.append([row[0],G,row[3],row[1],M,row[3]+'-'+M+'-Z2',S_LOW,len(S_LOW),R1])
        else: 
            
            OUTPUT1.append([row[0],G,row[3],row[1],M,row[3]+'-'+M+'-Z1','NO ZOND',0,0])


for G in GENE_LIST:
    df2 = pd.read_sql_query('SELECT DISTINCT("SEQ"),"Tm","Gc","GENE","MRNA","TYPE","RATE","LEN" FROM "main"."table" WHERE "GENE"=' + str(G) + ' AND ("MRNA"="ALL" OR "TYPE"="SOLO_GENE") AND "RATE" >60 ORDER BY "RATE" DESC,"LEN" DESC', blast_db)
    block=df1[df1[2]==G] 
    if df2.shape[0]>0:
        df2=df2.sort_values(by='LEN',ascending=False)
        ZA=df2.iloc[0]
        S_LOW=ZA['SEQ']
        R1=ZA['RATE']            
        OUTPUT1.append([block[0].iloc[0],G,block[3].iloc[0],'--','ALL',block[3].iloc[0] + '-ALL-Z1',S_LOW,len(S_LOW),R1])
        if df2.shape[0]>1:
            df2=df2.sort_values(by='LEN',ascending=False)
            ZA=df2.iloc[1]
            S_LOW=ZA['SEQ']
            R1=ZA['RATE']            
            OUTPUT1.append([block[0].iloc[0],G,block[3].iloc[0],'--','ALL',block[3].iloc[0] + '-ALL-Z2',S_LOW,len(S_LOW),R1])
        else:
            OUTPUT1.append([block[0].iloc[0],G,block[3].iloc[0],'--','ALL',block[3].iloc[0] + '-ALL-Z2',S_LOW,len(S_LOW),R1])
    else: 
            
        OUTPUT1.append([block[0].iloc[0],G,block[3].iloc[0],'--','ALL',block[3].iloc[0] + '-ALL-Z1','NO ZOND',0,0])                              


GENE_LIST=np.unique(df1a[6])
for G in GENE_LIST:
    block=df1a[df1a[6]==G]

    for index,row in block.iterrows():
        M=row[3]

        
        df2 = pd.read_sql_query('SELECT DISTINCT("SEQ"),"Tm","Gc","GENE","MRNA","TYPE","RATE","LEN","POS" from "main"."table" WHERE "GENE"="' + G + '" AND "MRNA"="' + str(M) + '"  AND "RATE" >60 ORDER BY "RATE" DESC,"LEN" DESC', blast_db)
        
        if df2.shape[0]>0:
            df2=df2.sort_values(by=['POS','LEN'],ascending=[False,False])
            ZA=df2.iloc[0]
            S_LOW=ZA['SEQ']
            R1=ZA['RATE']       
            OUTPUT1.append([row[0],row[3],'VIRUS',row[1],M,'VIRUS'+'-'+M+'-Z1',S_LOW,len(S_LOW),R1])
            if df2.shape[0]>1:
                df2=df2.sort_values(by=['POS','LEN'],ascending=[True,False])
                ZA=df2.iloc[0]
                S_LOW=ZA['SEQ']
                R1=ZA['RATE']       
                OUTPUT1.append([row[0],row[3],'VIRUS',row[1],M,'VIRUS'+'-'+M+'-Z1',S_LOW,len(S_LOW),R1])
            else:
                OUTPUT1.append([row[0],row[3],'VIRUS',row[1],M,'VIRUS'+'-'+M+'-Z2',S_LOW,len(S_LOW),R1])
        else: 
            
            OUTPUT1.append([row[0],row[3],'VIRUS',row[1],M,'VIRUS'+'-'+M+'-Z1','NO ZOND',0,0])

OUTPUT1=pd.DataFrame(OUTPUT1)
OUTPUT1.columns=['C1','GENE','GENE_NAME','C2','MRNA','LABEL','Z','L','R']
OUTPUT1=OUTPUT1.sort_values(by=['C1','C2'],ascending=[True,True])
OUTPUT1.to_csv('00_REZULT.csv',index=False,header=True,sep=';')


GOOD=OUTPUT1[OUTPUT1['Z']!='NO ZOND']
BAD=OUTPUT1[OUTPUT1['Z']=='NO ZOND']
NO_GENE=BAD[BAD['MRNA']=='ALL']
GENE=GOOD[GOOD['MRNA']=='ALL']
ISO=GOOD[GOOD['C2']!='--']
ALL=OUTPUT1

writer = pd.ExcelWriter('00_REZULT.xlsx')
ALL.to_excel(writer,index=False,sheet_name='ALL')
GOOD.to_excel(writer,index=False,sheet_name='GOOD')
ISO.to_excel(writer,index=False,sheet_name='ISO')
GENE.to_excel(writer,index=False,sheet_name='GENE')
BAD.to_excel(writer,index=False,sheet_name='BAD')
NO_GENE.to_excel(writer,index=False,sheet_name='NO_GENE')
writer.save() 