import pandas as pd
import numpy as np

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

import re

def QA(src):
    src=src.split(':')
    return src


def QA2(Z):
    try:
        blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02_VIRUS.db') 
        zond_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_01_VIRUS.db') 
        df1 = pd.read_sql_query('SELECT * FROM "main"."table" where qseqid=' + str(Z), blast_db) 
        S=df1['qseq'].iloc[0]
        L=len(S)
        df2 = pd.read_sql_query('SELECT * FROM "main"."table_f" where seq="' + S +'"', zond_db)         
        dfZ=df1[(df1['gaps']==0) & (df1['mismatch']==0) & (df1['length']>=0.7*L)]
  
        df2_Tm=np.unique(df2['Tm'])[0]
        df2_Gc=np.unique(df2['Gc'])[0]
    #        df2_Sl=np.unique(df2['sim_level'])[0]
    
    
        df2_Sl=0
        DST=[S,L,'--','--','--','--',df2_Tm,df2_Gc,df2_Sl,0,100]        
        if dfZ.shape[0]==0:
            DST=[S,L,'--','--','--','NO HITS TOTAL',df2_Tm,df2_Gc,df2_Sl,0,100] 
        else:
            l=pd.DataFrame(list(map(QA,list(dfZ['sseqid']))))
            GENE=list(np.unique(l[0]))
            GENE_S=[x for x in GENE if x != '-'] # num of target genes                               
        return DST        
    except:
        print('error on + ' + str(Z))


xlsx = pd.ExcelFile('SRC_DATA.xlsx')
REF_LIST = xlsx.parse('VIRUS')

blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02_VIRUS.db') 
df = pd.read_sql_query('SELECT DISTINCT(qseqid) FROM "main"."table"', blast_db)
df=list(df['qseqid'].values)
df=df[1:10]


Z=202689
start=time.time()
pool = ThreadPool(24)
DD=pool.map(QA2,df)
pool.close()
pool.join()
DD=pd.DataFrame(DD)
DD.columns=['SEQ','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
DD.to_csv('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.csv',index=False)
print(str(time.time()-start))


