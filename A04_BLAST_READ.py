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
    dst=[]
    src=src.split('|')
    ref=src[3]
    ref=ref.split('.')
    dst=ref[0]
    return dst


def QA2(Z):
    try:
        blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02.db') 
        zond_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_01.db')        
    
        query_b='SELECT seq,Tm,Gc,sim_level FROM "main"."table_f" where IDX=' + str(Z)        
        df2 = pd.read_sql_query(query_b, zond_db)        
        df2_Tm=np.unique(df2['Tm'])[0]
        df2_Gc=np.unique(df2['Gc'])[0]
        df2_Sl=np.unique(df2['sim_level'])[0]    
        L=len(np.unique(df2['seq'])[0])
        S=np.unique(df2['seq'])[0]
        S1='---'
        query_a='SELECT * FROM "main"."table" where qseqid=' + str(Z) + ' and gaps=0 and mismatch=0 and length>=' + str(0.7*L)
        dfZ = pd.read_sql_query(query_a,blast_db)          

        DST=[S,S1,L,'--','--','--','--',df2_Tm,df2_Gc,df2_Sl,0,100]        
        if dfZ.shape[0]==0:
            DST=[S,S1,L,'--','--','--','NO HITS TOTAL',df2_Tm,df2_Gc,df2_Sl,0,100] 
        else:
            hits=list(dfZ['sseqid'].values)
            hits_seq=list(dfZ['qseq'].values)
            hits=map(QA,hits)
            
            gen_list=[]
            mrna_list=[]
            h_list_a=[]
            h_list_b=[]
            S1x=[]
            for H,SH in zip(hits,hits_seq):
                S1x.append(SH) # what is exactly was aligned
                try:
                    pos=ML.index(H) # is in list of approved
                    gen_list.append(GL[pos])  # list of interscope gene 
                    mrna_list.append(H)  # list of interscope iso
                    h_list_a.append(H)
                except:
                    h_list_a.append(H) # list all iso
                    h_list_b.append(H) # list of outscope iso
            S1=str(np.unique(S1x))     
            if len(h_list_b)>0: #even 1 hit outside known isoform list
                m = [string for string in h_list_b if re.match("N(M|R)(_)", string)]   # list of refseq outlayers
                if len(m)>0: # there is refseq outlayers                
                    DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),str(m),'OUTSCOPE REFSEQ',df2_Tm,df2_Gc,df2_Sl,0,100] 
                else: #all outlayers are modeled
                    if len(np.unique(gen_list))==1: #only one gene
                       D=REF_LIST[REF_LIST[2]==np.unique(gen_list)[0]]
                       C=D.shape[0]
                       m0 = [string for string in mrna_list if re.match("N(M|R)(_)", string)] # list of refseq interscope
                       m1 = [string for string in h_list_b if re.match("X(M|R)(_)", string)] # list of model outscope
                       m2 = [string for string in mrna_list if re.match("X(M|R)(_)", string)] # list of refseq interscope
                       if len(m0)==1:
                           if len(m2)>0:# only one isoform              
                               DST=[S,S1,L,str(np.unique(gen_list)),str(m0),str(m1 + m2),'SOLO WITH INTERSCOPE AND OUTSCOPE MODEL',df2_Tm,df2_Gc,df2_Sl,1,C] 
                           else:
                               DST=[S,S1,L,str(np.unique(gen_list)),str(m0),str(m1),'SOLO WITH OUTSCOPE MODEL',df2_Tm,df2_Gc,df2_Sl,1,C] 
                       else:
                           if len(m2)>0:# only one isoform              
                               DST=[S,S1,L,str(np.unique(gen_list)),str(m0),str(m1 + m2),'GROUP WITH INTERSCOPE AND OUTSCOPE MODEL',df2_Tm,df2_Gc,df2_Sl,len(m0),C] 
                           else:
                               DST=[S,S1,L,str(np.unique(gen_list)),str(m0),str(m1),'GROUP WITH OUTSCOPE MODEL',df2_Tm,df2_Gc,df2_Sl,len(m0),C]  
                    else:
                        DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),'--','OUTSCOPE MULTIGENE',df2_Tm,df2_Gc,df2_Sl,0,100]   
            else:        #all hits are in scope  
                if len(gen_list)==0: #dump situation
                   DST=[S,S1,L,'--','--','--','NO HITS TOTAL',df2_Tm,df2_Gc,df2_Sl,0,100] 
                else:
                   if len(np.unique(gen_list))==1: #only one gene
                       D=REF_LIST[REF_LIST[2]==np.unique(gen_list)[0]]
                       C=D.shape[0]
                       if len(np.unique(mrna_list))==1: # only one isoform
                           DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),'--','SOLO',df2_Tm,df2_Gc,df2_Sl,1,C] 
                       else:
                           m = [string for string in mrna_list if re.match("N(M|R)(_)", string)]
                           m1 = [string for string in mrna_list if re.match("X(M|R)(_)", string)]
                           if len(m)==1: # list of isoform, but only one is N*
                               DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),str(m1),'SOLO WITH MODEL',df2_Tm,df2_Gc,df2_Sl,1,C]
                           else:
                               # more, then 1 isoform and more then 1 is N*
                               DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),'--','GROUP',df2_Tm,df2_Gc,df2_Sl,len(m),C] 
                   else: # more the 1 gene from scope in hits
                       DST=[S,S1,L,str(np.unique(gen_list)),str(np.unique(mrna_list)),'--','INTRASCOPE MULTIGENE',df2_Tm,df2_Gc,df2_Sl,0,100]
    
        return DST        
    except Exception as e:
        print('> error on  ' + str(Z) + ' | ' + 'query_a: ' + query_a + ' | ' + ' query b: ' + query_b)

xlsx = pd.ExcelFile('SRC_DATA.xlsx')
li=xlsx.sheet_names
i=li[0]
df = xlsx.parse('HUMAN')
REF_LIST=df


 
DATA_2 = pd.read_csv('REG_LIST.csv')
GL=np.array(DATA_2['0'].values)
ML=list(DATA_2['1'].values)
#MLA=np.array(DATA_2['1'].values)

blast_db = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02.db') 
df = pd.read_sql_query('SELECT DISTINCT(qseqid)  FROM "main"."table"', blast_db)

id_list=list(df['qseqid'].values)


start=time.time()
pool = ThreadPool(24)
DD=pool.map(QA2,id_list)
pool.close()
pool.join()
DD=pd.DataFrame(DD)
DD.columns=['SEQ','SEQ1','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
DD.to_csv('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03.csv',index=False)
print(str(time.time()-start))


csv_database = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_03.db')       
dfa=['aaa','aaa','aaa','aaa','aaa','aaa','aaa',1.1999,1.1999,1.199,0,0]
dfa=pd.DataFrame(dfa)
dfa=dfa.transpose()
dfa.columns=['SEQ','SEQ1','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
dfa.to_sql('table', csv_database, if_exists='replace') 


chunksize = 16384*8
for df in pd.read_csv('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03.csv', chunksize=chunksize, iterator=True,sep=',',low_memory=False):         
      df.columns=['SEQ','SEQ1','LEN','HIT_GENE','HIT_MRNA','OUTSCOPE','TYPE','Tm','Gc','SIM_LEVEL','H1','H2']
      df = df.rename(columns={c: c.replace(' ', '') for c in df.columns}) 
      df.to_sql('table', csv_database, if_exists='append')  

              
              
db1 = sqlite3.connect('/home/www_adm/ssd/CHIP_OUTPUT/STEP_03.db')
cursor = db1.cursor()
try:
    sql = 'CREATE INDEX ix1 ON "table" ("SEQ")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix1a ON "table" ("SEQ1")'
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