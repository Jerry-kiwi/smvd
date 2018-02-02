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

import os   
TAG='/home/www_adm/ssd/CHIP_OUTPUT/BLAST_FASTA'
shutil.rmtree(TAG, ignore_errors=True)
os.mkdir(TAG)
    
csv_database = create_engine('sqlite:////home/www_adm/ssd/CHIP_OUTPUT/STEP_01.db')       
df = pd.read_sql_query('SELECT DISTINCT(seq),IDX FROM "main"."table_f" order by IDX', csv_database) 

seq_L=list(df['seq'])
seq_L=np.array(seq_L)

id_L=list(df['IDX'])
id_L=np.array(id_L)

CHUNK_SIZE=1000
SX=[]
start=0;
finish=start+CHUNK_SIZE
while finish<=len(seq_L)+CHUNK_SIZE-2:
    SX.append([start,finish])
    start=finish+1
    finish=start+CHUNK_SIZE
SX=np.array(SX)
c=0
for CK in range(0,len(SX)):
    seq=seq_L[SX[CK,0]:SX[CK,1]+1]
    idx=id_L[SX[CK,0]:SX[CK,1]+1]
    rec_list=[]
    for s,ids in zip(seq,idx):
        rec = SeqRecord(Seq(s, generic_nucleotide),id=str(ids),description='')
        rec_list.append(rec)    
        c=c+1
    SeqIO.write(rec_list, '/home/www_adm/ssd/CHIP_OUTPUT/BLAST_FASTA/' + str(CK) + '.fasta', "fasta")            


