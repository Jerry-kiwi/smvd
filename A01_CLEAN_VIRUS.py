import pandas as pd
import numpy as np

from Bio import Entrez,SeqUtils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_nucleotide
from Bio.SeqUtils import MeltingTemp

import itertools
from itertools import count # izip for maximum efficiency

import ahocorasick
import re

import multiprocessing
from multiprocessing import Pool as ThreadPool
import time


import os
import subprocess
import shutil

from sqlalchemy import create_engine
import sqlite3

MIN_ZOND_LEN=20
MAX_ZOND_LEN=30
EXON_CMP=.66

MIN_TM=40
MAX_TM=70
MIN_GC=40
MAX_GC=70

MAX_REPEATS=4

def get_overlapped_chunks(textin, chunksize, overlapsize):  
    return [ textin[a:a+chunksize] for a in range(0,len(textin), chunksize-overlapsize)]

def check_phys(src):  
    
    dst=0
    TM=MeltingTemp.Tm_NN(src)
    GC=SeqUtils.GC(src)
    regex = r"(C|G)\1{3,}"
    HOMO=re.findall(regex, src) # 3 or more
    
    regex = r"(A|T)\1{4,}"
    HOMO2=re.findall(regex, src) # 3 or more
    
    if ((TM<MIN_TM) | (TM>MAX_TM)): dst=1
    if ((GC<MIN_GC) | (GC>MAX_GC)): dst=2
    if (len(HOMO)>0): dst=3
    if (len(HOMO2)>0): dst=4
    response=[TM,GC,len(HOMO),dst]
    return response
    
def acho(src):
    dst=[src,src,len(src)/len(src)]
#    tmp=[]
#    for j in range(10,len(src)+1):
#        a=src[0:j]
#        b=src[len(src)-j:len(src)]
#        if (len(a)>0):tmp.append(a)
#        if (len(b)>0):tmp.append(b)   
#
#    A=ahocorasick.Automaton() 
#    c=0
#    for i in tmp:
#       A.add_word(i, (c,i))
#       c=c+1
#       
#    A.make_automaton()
#    oc=[]
#    seq=[]
#    for item in A.iter(GENE_SEQ):
#        oc.append(len(item[1][1]))
#        seq.append(item[1][1])
#    oc=np.array(oc)/len(src)
#    try:
#        dst=[src,seq[np.argmax(oc)],oc[np.argmax(oc)]]
#    except:
#        dst=['','',0]

    return dst


csv_database_A = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_01_VIRUS.db')       
dfa=['XXX',20,1,1,1,1,1]
dfa=pd.DataFrame(dfa)
dfa=dfa.transpose()
dfa.columns=['seq','L','sim_level','Tm','Gc','dG','sp']
dfa.to_sql('table', csv_database_A, if_exists='replace')
        
        
        
xlsx = pd.ExcelFile('SRC_DATA.xlsx')
df = xlsx.parse('VIRUS')

ID_A=np.unique(df[0])
for I in ID_A:
    GENE_BLOCK=df[df[0]==I]
    MRNA_COUNT=GENE_BLOCK.shape[0]
    for M in range(1,MRNA_COUNT+1):
        MRNA_BLOCK=GENE_BLOCK[GENE_BLOCK[1]==M]

        start = time.time()
        Gp=MRNA_BLOCK[0].iloc[0]
        Mp=MRNA_BLOCK[1].iloc[0]
        
        GENE_FILE_NAME='/home/www_adm/ownCloud/NNIIEM_SCIENCE/CELL_DEATH/A04_CHIP_REV_2018/GENE_FASTA/' + str(Gp) + '_0.fasta'
        MRNA_FILE_NAME='/home/www_adm/ownCloud/NNIIEM_SCIENCE/CELL_DEATH/A04_CHIP_REV_2018/MRNA_FASTA/'+ str(Gp) + '_' + str(Mp) + '.fasta'
        
        for record in SeqIO.parse(GENE_FILE_NAME, "fasta"):          
            GENE_SEQ=str(record.seq)
        
        for record in SeqIO.parse(MRNA_FILE_NAME, "fasta"):   
            MRNA_SEQ=str(record.seq)
        
        MRNA_SEQ=MRNA_SEQ.upper()

#            
        VOC_LIST=[]
        for ZOND_LEN in range(MIN_ZOND_LEN,MAX_ZOND_LEN+1):
            chunks = get_overlapped_chunks(MRNA_SEQ, ZOND_LEN, ZOND_LEN-1)
            VOC_LIST=VOC_LIST+chunks
        
        
        VOC_LIST=list(set(VOC_LIST))
        pos_good=np.array([i for i, j in zip(count(), VOC_LIST) if len(j) >=MIN_ZOND_LEN])               
        VOC_LIST = [VOC_LIST[index] for index in pos_good] # зонды с длиной более 20 н.о. (особенности алгоритма нарезания)

        rec_list=[]
        for V in VOC_LIST:
            rec = SeqRecord(Seq(V, generic_nucleotide),id=MRNA_BLOCK[5].iloc[0],description='')
            rec_list.append(rec)
        SeqIO.write(rec_list, '/home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' + str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.fasta', "fasta")            
        
        cmd='hybrid-ss-min /home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' + str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.fasta' + ' -E -t 42 -T 42.1 -n DNA -o /home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' +  str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0])
        subprocess.check_call(cmd,shell=True)
        dG=pd.read_csv('/home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' +  str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.dG', sep='\t')
        dG.columns=['a','b','c']
        
        os.remove('/home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' +  str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.dG')
        os.remove('/home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' +  str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.run')
        os.remove('/home/www_adm/ssd/CHIP_OUTPUT/DRAFT/' +  str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + '.fasta')

        
        pool = ThreadPool(multiprocessing.cpu_count())
        V1=pool.map(check_phys,VOC_LIST)
        pool.close()
        pool.join()
        
        V1X=pd.DataFrame(V1)
        V1Y=pd.concat([pd.DataFrame(VOC_LIST),V1X,dG['b']],axis=1) #block with  phys and dG
        
       
        V1Y.columns=['S','T','G','H','F','E']

        posA=np.array(V1Y[(V1Y['F']==0) & (V1Y['E']>=-3) & (V1Y['E']<=10)].index)
        
        V1Y=V1Y.iloc[posA,:]
        VOC_LIST=list(V1Y['S'])
        
        L = pd.DataFrame(list(map(len, VOC_LIST)))

        pool = ThreadPool(multiprocessing.cpu_count())
        V1A=pool.map(acho,VOC_LIST)
        pool.close()
        pool.join()          
        
        V1A=pd.DataFrame(V1A)
        V1A.columns=['zond','zond_part','sim_to_gene']
        V1A=V1A[['zond_part','sim_to_gene']]
        
        V1Y=V1Y.reset_index(drop=True)
        V1A=V1A.reset_index(drop=True)
        L=L.reset_index(drop=True)
        RESPONSE=pd.concat([V1Y,L,V1A],axis=1)
        RESPONSE.columns=['seq','Tm','Gc','H','F','dG','L','sp','sim_level']
        RESPONSE=RESPONSE[['seq','L','sim_level','Tm','Gc','dG','sp']]
        RESPONSE.to_sql('table', csv_database_A, if_exists='append')
        
        print(str(MRNA_BLOCK[0].iloc[0]) + '_' + str(MRNA_BLOCK[1].iloc[0]) + ' ' + MRNA_BLOCK[5].iloc[0] + ' ' + str(len(VOC_LIST)) + ' ' + str(time.time() - start))


#CREATE TABLE `table_2` (
#	`IDX`	INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE,
#	`seq`	TEXT,
#	`gene`	INTEGER,
#	`L`	INTEGER,
#	`Tm`	REAL,
#	`Gc`	REAL,
#	`sim_level`	REAL
#);
#create table "table_f" as select distinct(seq),gene,L,Tm,Gc,sim_level from "table"     
#INSERT INTO "table_2" SELECT NULL,seq,gene,L,Tm,Gc,sim_level FROM "table_f"
#CREATE INDEX "ix1" ON "table_2" ("seq")
#CREATE INDEX "ix2" ON "table_2" ("IDX")
#CREATE INDEX "ix3" ON "table_2" ("seq","IDX")
#DROP TABLE "table_f"
#ALTER TABLE "table_2" RENAME TO "table_f"        

#sys.stdout = orig_stdout
#f.close()
        

