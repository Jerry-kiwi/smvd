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
import re



#import sys
#orig_stdout = sys.stdout
#f = open('000.txt', 'w')
#sys.stdout = f


def sel_human(input_param):
    SOLO_BLOCK=[]
    G=input_param[0]
    M=input_param[1]

    df2 = pd.read_sql_query('SELECT * from "main"."table" where SIM_LEVEL<=.7 and HIT_GENE="[' + str(G) + ']" and HIT_MRNA LIKE "%' + str(M) + '%"' , blast_db)     
    SOLO_BLOCK=df2[(df2['TYPE']=="SOLO") | (df2['TYPE']=="SOLO WITH MODEL") | (df2['TYPE']=="SOLO WITH OUTSCOPE MODEL") | (df2['TYPE']=='SOLO WITH INTERSCOPE AND OUTSCOPE MODEL')]
    SOLO_BLOCK=SOLO_BLOCK[['SEQ','Tm','Gc']]
    SOLO_BLOCK['GENE']=G
    SOLO_BLOCK['MRNA']=M
    SOLO_BLOCK['TYPE']='SOLO'

    if SOLO_BLOCK.shape[0]>0:
        SOLO_BLOCK=SOLO_BLOCK
    else:
        dfa=['aaa','aaa','aaa','aaa','DUMP','aaa']
        dfa=pd.DataFrame(dfa)
        dfa=dfa.transpose()
        dfa.columns=['SEQ','Tm','Gc','GENE','MRNA','TYPE']
        SOLO_BLOCK=dfa

    SOLO_G_BLOCK=[]

    df2 = pd.read_sql_query('SELECT * from "main"."table" where SIM_LEVEL=1 and H2=1 and H1=1 and HIT_GENE="[' + str(G) + ']" and HIT_MRNA LIKE "%' + str(M) + '%"' , blast_db)     
    SOLO_G_BLOCK=df2[(df2['TYPE']=="SOLO") | (df2['TYPE']=="SOLO WITH MODEL") | (df2['TYPE']=="SOLO WITH OUTSCOPE MODEL") | (df2['TYPE']=='SOLO WITH INTERSCOPE AND OUTSCOPE MODEL')]
    SOLO_G_BLOCK=SOLO_G_BLOCK[['SEQ','Tm','Gc']]
    SOLO_G_BLOCK['GENE']=G
    SOLO_G_BLOCK['MRNA']=M
    SOLO_G_BLOCK['TYPE']='SOLO_GENE'

    if SOLO_G_BLOCK.shape[0]>0:
        SOLO_G_BLOCK=SOLO_G_BLOCK
    else:
        dfa=['aaa','aaa','aaa','aaa','DUMP','aaa']
        dfa=pd.DataFrame(dfa)
        dfa=dfa.transpose()
        dfa.columns=['SEQ','Tm','Gc','GENE','MRNA','TYPE']
        SOLO_G_BLOCK=dfa

    GROUP_BLOCK=[]

    df2 = pd.read_sql_query('SELECT * from "main"."table" where HIT_GENE="[' + str(G) + ']"' , blast_db)  
    GROUP_BLOCK=df2[(df2['TYPE']=="GROUP") | (df2['TYPE']=="GROUP WITH OUTSCOPE MODEL") | (df2['TYPE']=='GROUP WITH INTERSCOPE AND OUTSCOPE MODEL')]
    GROUP_BLOCK=GROUP_BLOCK[['SEQ','Tm','Gc','H1','H2']]
    
    GROUP_BLOCK=GROUP_BLOCK[GROUP_BLOCK['H1']>0.6*GROUP_BLOCK['H2']]
    
    GROUP_BLOCK['GENE']=G
    GROUP_BLOCK['MRNA']='ALL'
    GROUP_BLOCK['TYPE']='GENE'
    GROUP_BLOCK=GROUP_BLOCK[['SEQ','Tm','Gc','GENE','MRNA','TYPE']] 

    if GROUP_BLOCK.shape[0]>0:
        GROUP_BLOCK=GROUP_BLOCK

    else:
        dfa=['aaa','aaa','aaa','aaa','DUMP','aaa']
        dfa=pd.DataFrame(dfa)
        dfa=dfa.transpose()
        dfa.columns=['SEQ','Tm','Gc','GENE','MRNA','TYPE']
        GROUP_BLOCK=dfa
    return pd.concat([SOLO_BLOCK,SOLO_G_BLOCK,GROUP_BLOCK],axis=0)


    
def sel_virus(input_param):
    SOLO_BLOCK=[]
    G=input_param[0]
    M=input_param[1]
    df2 = pd.read_sql_query('SELECT * from "main"."table" where HIT_GENE="' + G + '" and HIT_MRNA LIKE "%' + M + '%"' , blast_db)    
    SOLO_BLOCK=df2[df2['TYPE']=="SOLO_VIRUS"]
    SOLO_BLOCK=SOLO_BLOCK[['SEQ','Tm','Gc']]
    SOLO_BLOCK['GENE']=G
    SOLO_BLOCK['MRNA']=M
    SOLO_BLOCK['TYPE']='SOLO'

    if SOLO_BLOCK.shape[0]>0:
        SOLO_BLOCK=SOLO_BLOCK.drop_duplicates(subset='SEQ')
        return SOLO_BLOCK
    else:
        dfa=['aaa','aaa','aaa','aaa','DUMP','aaa']
        dfa=pd.DataFrame(dfa)
        dfa=dfa.transpose()
        dfa.columns=['SEQ','Tm','Gc','GENE','MRNA','TYPE']
        return dfa


def gen_mrna(input_param):
    SOLO_BLOCK=[]
    GN=input_param[0]
    MN=input_param[1]
    for record in SeqIO.parse('/home/www_adm/ownCloud/NNIIEM_SCIENCE/CELL_DEATH/A04_CHIP_REV_2018/MRNA_FASTA/' + str(GN) + '_' + str(MN) + '.fasta', "fasta"):          
        SEQ=str(record.seq)
    SEQ=SEQ.upper()
    SOLO_BLOCK.append([GN,MN,SEQ])
    SOLO_BLOCK=pd.DataFrame(SOLO_BLOCK)
    SOLO_BLOCK.columns=['G','M','S']
    return SOLO_BLOCK


def cl_len(mrna_id,mrna_seq):
    try:
        if mrna_id!='ALL':
            pos=df1[df1[5]==mrna_id]
            GN=int(pos[0].values[0])
            MN=int(pos[1].values[0])
            SEQ=SEQ_BLOCK['S'][(SEQ_BLOCK['G']==GN) & (SEQ_BLOCK['M']==MN)]
            SEQ=SEQ[0]
            L=len(SEQ)
            pos_list=[m.start() for m in re.finditer(mrna_seq, SEQ)] 
            pos_list=np.array(pos_list)
            #if len(pos_list)>1:
            #    print ('many points ' + str(pos_list) + ' with ' + str(GN) + '_' + str(MN) + '  ' + mrna_id + '  ' + mrna_seq)
            pos_list=pos_list[0]*100/L
            pos_list=np.ceil(pos_list)
        else:
            pos_list=0
        return pos_list
    except:
        try:
            if mrna_id!='ALL':
                pos=df2[df2[3]==mrna_id]
                GN=int(pos[0].values[0])
                MN=int(pos[1].values[0])
                SEQ=SEQ_BLOCK['S'][(SEQ_BLOCK['G']==GN) & (SEQ_BLOCK['M']==MN)]
                SEQ=SEQ[0]
                L=len(SEQ)
                pos_list=[m.start() for m in re.finditer(mrna_seq, SEQ)] 
                pos_list=np.array(pos_list)
                #if len(pos_list)>1:
                #    print ('many points ' + str(pos_list) + ' with ' + str(GN) + '_' + str(MN) + '  ' + mrna_id + '  ' + mrna_seq)
                pos_list=pos_list[0]*100/L
                pos_list=np.ceil(pos_list)
            
            else:
                pos_list=0
            return pos_list
        except Exception as e:
            print ('error with ' + mrna_id + '  ' + mrna_seq + '   ' + str(e))

start=time.time()
print('Reading src human xlxs')
xlsx = pd.ExcelFile('SRC_DATA.xlsx')
li=xlsx.sheet_names
df1 = xlsx.parse('HUMAN')
blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_03.db') 
input_param=[]
glist=np.unique(df1[2])
for G in glist:
    block=df1[df1[2]==G]
    for index2,row2 in block.iterrows():
        input_param.append([row2[2],row2[5]])       
print(str(time.time()-start))  
print('Finish src human preparing')

         
start=time.time()
print('Gen src human block')
pool = ThreadPool(24)
OUTPUTA=pool.map(sel_human,input_param)
pool.close()
pool.join()
print(str(time.time()-start))  
print('Finish src human solo')


start=time.time()
print('Reading src virus xlxs')
xlsx = pd.ExcelFile('SRC_DATA.xlsx')
li=xlsx.sheet_names
df1 = xlsx.parse('VIRUS')

blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.db') 


input_param=[]
glist=np.unique(df1[2])
for G in glist:
    block=df1[df1[2]==G]
    for index2,row2 in block.iterrows():
        input_param.append([row2[6],row2[3]])      

print(str(time.time()-start))  
print('Finish solo virus preparing')

start=time.time()
print('Gen solo virus block')
pool = ThreadPool(24)
OUTPUTC=pool.map(sel_virus,input_param)
pool.close()
pool.join()
print(str(time.time()-start))  
print('Finish src virus')


OUTPUTC=pd.concat(OUTPUTC)       
#OUTPUTC=OUTPUTC.drop_duplicates(subset='SEQ')
#
OUTPUTA=pd.concat(OUTPUTA)
#OUTPUTA=OUTPUTA.drop_duplicates(subset='SEQ')

OUTPUT1=pd.concat([OUTPUTA,OUTPUTC])
OUTPUT1=OUTPUT1[OUTPUT1['MRNA']!='DUMP']

start=time.time()
print('Gen rate')
TM_V=np.array(list(OUTPUT1['Tm']))
GC_V=np.array(list(OUTPUT1['Gc']))


x=np.arange(min(GC_V),max(GC_V),0.1)
mu=np.mean(x)
sig=5
y=(1/(sig*math.sqrt(2*math.pi)))*np.exp(-1*(np.power((x-mu),2)/(2*math.pow(sig,2))))
GC_RATE=(y-min(y))/(max(y)-min(y));
GC_RATE=np.interp(GC_V, x,GC_RATE)

x=np.arange(min(TM_V),max(TM_V),0.1)
mu=np.mean(x)
sig=3
y=(1/(sig*math.sqrt(2*math.pi)))*np.exp(-1*(np.power((x-mu),2)/(2*math.pow(sig,2))))
TM_RATE=(y-min(y))/(max(y)-min(y));
TM_RATE=np.interp(TM_V, x,TM_RATE)


FINAL_RATE=TM_RATE*0.7+GC_RATE*0.3
FINAL_RATE=FINAL_RATE*100
FINAL_RATE=np.floor(FINAL_RATE)
OUTPUT1['RATE']=FINAL_RATE
OUTPUT1['LEN']=list(map(len,OUTPUT1['SEQ']))
print(str(time.time()-start))  
print('Finish rate')


start=time.time()
print('Reading src for pos block xlxs')
xlsx = pd.ExcelFile('SRC_DATA.xlsx')
df1 = xlsx.parse('HUMAN')
df2 = xlsx.parse('VIRUS')


input_param1=[]
glist=np.unique(df1[2])
for G in glist:
    block=df1[df1[2]==G]
    for index2,row2 in block.iterrows():
        input_param1.append([int(row2[0]),int(row2[1])])      

input_param2=[]
glist=np.unique(df2[2])
for G in glist:
    block=df2[df2[2]==G]
    for index2,row2 in block.iterrows():
        input_param2.append([row2[0],row2[1]])      

       
input_param=input_param1+input_param2

start=time.time()
print('Gen seq block')
pool = ThreadPool(24)
CL_LIST=pool.map(gen_mrna,input_param)
pool.close()
pool.join() 
SEQ_BLOCK=pd.concat(CL_LIST)
print(str(time.time()-start))  
print('Finish seq block')

start=time.time()
print('Gen pos vec')
pool = ThreadPool(24)
CL_LIST=pool.starmap(cl_len,zip(OUTPUT1['MRNA'],OUTPUT1['SEQ']))
pool.close()
pool.join() 
OUTPUT1['POS']=CL_LIST
print(str(time.time()-start))  
print('Finish len calc')

        
start=time.time()
print('Gen db')       
csv_database = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_04b.db')       
OUTPUT1.to_sql('table', csv_database, if_exists='replace')  
print(str(time.time()-start))  
print('Finish db')
              
db1 = sqlite3.connect('/home/www_adm/ssd/CHIP_OUTPUT/STEP_04b.db')
cursor = db1.cursor()
try:
    sql = 'CREATE INDEX ix1 ON "table" ("SEQ")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix4 ON "table" ("GENE")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix5 ON "table" ("MRNA")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix7 ON "table" ("TYPE")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix8 ON "table" ("Tm")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix9 ON "table" ("Gc")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix10 ON "table" ("RATE")'
    cursor.execute(sql)

except:
    g=0
db1.close()
           