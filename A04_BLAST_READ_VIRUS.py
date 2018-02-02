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
        query_b='SELECT seq,Tm,Gc FROM "main"."table_f" where IDX=' + str(Z)
        
        df2 = pd.read_sql_query(query_b, zond_db)        
        df2_Tm=np.unique(df2['Tm'])[0]
        df2_Gc=np.unique(df2['Gc'])[0]
        df2_Sl=0  
        L=len(np.unique(df2['seq'])[0])
        S=np.unique(df2['seq'])[0]
        query_a='SELECT * FROM "main"."table" where qseqid=' + str(Z) + ' and gaps=0 and mismatch=0 and length>=' + str(0.7*L)
        dfZ = pd.read_sql_query(query_a,blast_db)          

        DST=[S,L,'--','--','--','--',df2_Tm,df2_Gc,df2_Sl,0,100]        
        if dfZ.shape[0]==0:
            DST=[S,L,'--','--','--','NO HITS TOTAL',df2_Tm,df2_Gc,df2_Sl,0,100] 
        else:
            l=pd.DataFrame(list(map(QA,list(dfZ['sseqid']))))
            GENE=list(np.unique(l[0]))
            GENE_S=[x for x in GENE if x != '-'] # num on target genes
            if len(GENE_S)>1:
               DST=[S,L,str(GENE),'--','--','INTERSCOPE MULTIGENE',df2_Tm,df2_Gc,df2_Sl,0,100]  
            else:
               REFERENCE_GENE=list(REF_LIST[6][REF_LIST[3]==GENE_S[0]].values)[0]
               outscope=l[l[1]!=REFERENCE_GENE]
               
               if outscope.shape[0]>0:
                   DST=[S,L,REFERENCE_GENE,'--',str(list(np.unique(outscope[1]))),'OUTSCOPE MULTIGENE',df2_Tm,df2_Gc,df2_Sl,0,100] 
               else:
                   marker='-:'+REFERENCE_GENE+':0:0'
                   DX=dfZ[dfZ['sseqid']==marker] # hits from genome
                   DYa=l[l[0]==GENE_S[0]]
                   low=int(DYa[2].values[0])
                   high=int(DYa[3].values[0])
                   hit_list=[]
                   for index,row in DX.iterrows():
                       st=row['sstart']
                       fn=row['send']
                       la=min([st,fn])
                       lb=max([st,fn])
                       if la>low:
                           if lb<high:
                               hit_list.append(0)
                           else:
                               hit_list.append(0)
                       else:
                           hit_list.append(0)
                   hit_list=np.array(hit_list)
                   if sum(hit_list)==0:
                       DST=[S,L,REFERENCE_GENE,GENE_S[0],'--','SOLO_VIRUS',df2_Tm,df2_Gc,df2_Sl,0,100]
                   else:
                       DST=[S,L,REFERENCE_GENE,GENE_S[0],'--','SOLO_VIRUS WITH OUTGENE LOCATION',df2_Tm,df2_Gc,df2_Sl,0,100]
                               
        return DST        
    except Exception as e:
        print('> error on  ' + str(Z) + ' | ' + 'query_a: ' + query_a + ' | ' + ' query b: ' + query_b)


xlsx = pd.ExcelFile('SRC_DATA.xlsx')
REF_LIST = xlsx.parse('VIRUS')
 

blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02_VIRUS.db') 
df = pd.read_sql_query('SELECT DISTINCT(qseqid) FROM "main"."table"', blast_db)

id_list=list(df['qseqid'].values)


start=time.time()
pool = ThreadPool(24)
DD=pool.map(QA2,id_list)
pool.close()
pool.join()
DD=pd.DataFrame(DD)
DD.columns=['SEQ','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
DD.to_csv('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.csv',index=False)
print(str(time.time()-start))


csv_database = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.db')       
dfa=['aaa','aaa','aaa','aaa','aaa','aaa',1.1999,1.1999,1.199,0,0]
dfa=pd.DataFrame(dfa)
dfa=dfa.transpose()
dfa.columns=['SEQ','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
dfa.to_sql('table', csv_database, if_exists='replace') 


chunksize = 16384*4
for df in pd.read_csv('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.csv', chunksize=chunksize, iterator=True,sep=',',low_memory=False):         
      df.columns=['SEQ','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
      df = df.rename(columns={c: c.replace(' ', '') for c in df.columns}) 
      rerun=True        
      while (rerun==True):
          try:
              df.to_sql('table', csv_database, if_exists='append')  
              rerun=False
          except Exception as e:
              time.sleep(3)
              
              
db1 = sqlite3.connect('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03_VIRUS.db')
cursor = db1.cursor()
try:
    sql = 'CREATE INDEX ix1 ON "table" ("SEQ")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix4 ON "table" ("HIT_GENE")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix5 ON "table" ("HIT_MRNA")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix6 ON "table" ("OUTSCOPE")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix7 ON "table" ("TYPE")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix8 ON "table" ("Tm")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix9 ON "table" ("Gc")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix10 ON "table" ("SIM_LEVEL")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix13 ON "table" ("H1")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix14 ON "table" ("H2")'
    cursor.execute(sql)
except:
    g=0
db1.close()