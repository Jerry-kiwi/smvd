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



def gen_fast(param):
    G=param
    input_fasta='/home/www_adm/ssd/CHIP_OUTPUT/BLAST_FASTA_VIRUS/' + str(G) + '.fasta'
    input_base='/home/www_adm/ssd/BLAST/VIRUS_CUSTOM/VIRUS_DB'
    output_report='/home/www_adm/ssd/CHIP_OUTPUT/BLAST_REPORT/' + str(G) + '.report'
    cmd='blastn -task megablast -word_size 16 -ungapped -qcov_hsp_perc 75 -perc_identity 100 -db '+ input_base + ' -query '+ input_fasta + ' -out ' + output_report + ' -outfmt "6 qseqid qstart qend qseq sseqid sstart send length mismatch gaps pident nident" -num_threads 1'
    subprocess.check_call(cmd,shell=True)
    
    chunksize = 8192
    for df in pd.read_csv(output_report, chunksize=chunksize, iterator=True,header=None,sep='\t'):         
          df.columns=['qseqid','qstart','qend','qseq','sseqid','sstart','send','length','mismatch','gaps','pident','nident']
          df = df.rename(columns={c: c.replace(' ', '') for c in df.columns}) 
          rerun=True        
          while (rerun==True):
              try:
                  df.to_sql('table', csv_database, if_exists='append')  
                  rerun=False
              except Exception as e:
                  time.sleep(3)
    os.remove(output_report)


csv_database = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_02_VIRUS.db')       
dfa=['aaa',1,1,'sss','aaa',1,1,1,1,1,1,1]
dfa=pd.DataFrame(dfa)
dfa=dfa.transpose()
dfa.columns=['qseqid','qstart','qend','qseq','sseqid','sstart','send','length','mismatch','gaps','pident','nident']
dfa.to_sql('table', csv_database, if_exists='replace') 


#db1 = sqlite3.connect('/home/www_adm/ssd/CHIP_OUTPUT/STEP_02.db')
#cursor = db1.cursor()
#cursor.execute("PRAGMA journal_mode=WAL")
#db1.close()
shutil.rmtree('/home/www_adm/ssd/CHIP_OUTPUT/BLAST_REPORT', ignore_errors=True)
os.mkdir('/home/www_adm/ssd/CHIP_OUTPUT/BLAST_REPORT')   
number_files = len(os.listdir('/home/www_adm/ssd/CHIP_OUTPUT/BLAST_FASTA_VIRUS/'))
GX=list(np.arange(0,number_files))
pool = ThreadPool(24)
V1=pool.map(gen_fast,GX)
pool.close()
pool.join()



db1 = sqlite3.connect('/home/www_adm/ssd/CHIP_OUTPUT/STEP_02_VIRUS.db')
cursor = db1.cursor()
try:
    sql = 'DELETE from "table" where qseqid="aaa"'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix2 ON "table" ("qseqid","qseq")'
    cursor.execute(sql)
    sql = 'CREATE INDEX ix3 ON "table" ("qseqid","gaps","mismatch","length")'
    cursor.execute(sql)
except:
    g=0
db1.close()