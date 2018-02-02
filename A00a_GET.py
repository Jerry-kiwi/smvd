import pandas as pd
import numpy as np

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as generic_nucleotide

import time


import os   
import shutil
TAG='/home/www_adm/ownCloud/NNIIEM_SCIENCE/CELL_DEATH/A04_CHIP_REV_2018/GENE_FASTA'
shutil.rmtree(TAG, ignore_errors=True)
os.mkdir(TAG)

import os   
TAG='/home/www_adm/ownCloud/NNIIEM_SCIENCE/CELL_DEATH/A04_CHIP_REV_2018/MRNA_FASTA'
shutil.rmtree(TAG, ignore_errors=True)
os.mkdir(TAG)

xlsx = pd.ExcelFile('INIT.xlsx')
li=xlsx.sheet_names
i=li[0]
df = xlsx.parse(i)

g_i_S=list(df['ID'])
g_n_S=list(df['GENE_NAME'])
f_n_S=list(df['DESC'])

count_list_S=np.arange(1,len(g_n_S)+1,1)

MAIN_BLOCK=[]
CHUNK_SIZE=30
SX=[]
start=0;
finish=start+CHUNK_SIZE
while finish<=len(g_i_S)+CHUNK_SIZE-4:
    SX.append([start,finish])
    start=finish+1
    finish=start+CHUNK_SIZE

SX=np.array(SX) 
TOTAL_M=[]
TOTAL_G=[]
TOTAL_C=0

VOC_LIST=[]
for CK in range(0,len(SX)):
    g_i=g_i_S[SX[CK,0]:SX[CK,1]+1]
    g_n=g_n_S[SX[CK,0]:SX[CK,1]+1]
    f_n=f_n_S[SX[CK,0]:SX[CK,1]+1]
    count_list=count_list_S[SX[CK,0]:SX[CK,1]+1]
    print(count_list)
    
    dx=[]
    Entrez.email = "solntsev.l.a@gmail.com"
    
    BLOCK_MAIN=[]
    
    rerun=True        
    while (rerun==True):
        try:
            handle = Entrez.efetch(db="gene", id=str(g_i), rettype="docsum",retmode="xml")   
            rerun=False
            record = Entrez.read(handle)
        except Exception as e:
            rerun=True
            print('Error getting info on chunk ' + str(CK) )
            time.sleep(20)

    OUTPUT=[]

    record=record['DocumentSummarySet']['DocumentSummary']
    for R1,count1 in zip(record,range(0,len(record))):
        print('Gene: ' + g_n[count1] + ' ' + str(g_i[count1]))
        try:
            from_pos=int(R1['LocationHist'][0]['ChrStart'])+1
            to_pos=int(R1['LocationHist'][0]['ChrStop'])+1
            DNA=str(R1['LocationHist'][0]['ChrAccVer'])
        except Exception as e:
            DNA=str(R1['GenomicInfo'][0]['ChrAccVer'])
            from_pos=int(R1['GenomicInfo'][0]['ChrStart'])+1
            to_pos=int(R1['GenomicInfo'][0]['ChrStop'])+1
#        from_pos=int(R1['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from'])+1
#        to_pos=int(R1['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to'])+1
        
#        try:
#            strand_text=R1['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand']
#            strand_text=strand_text['Na-strand'].attributes['value']
#        except Exception:
#            strand_text='not mRNA'     
#        
#    
#        if (strand_text=='plus'):
#            strand=1
#        else:
#            strand=2
        if (from_pos>to_pos):
            t=from_pos
            from_pos=to_pos
            to_pos=t
            strand=2
            strand_text='minus'
        else:
            strand=1
            strand_text='plus' 
            
        from_pos=min(from_pos,to_pos)  
        to_pos=max(from_pos,to_pos)
        
        zzz_NC_NAME=str(DNA) # номер хромосомы
        print(zzz_NC_NAME + ' ' + str(from_pos) + ':' + str(to_pos))        
        zzz_NC_NAME=zzz_NC_NAME.split('.');
        zzz_NC_NAME=zzz_NC_NAME[0]
    
        # strand - what strand of DNA to show (1 = plus or 2 = minus)
        # Для strand=2 номера координат в Gene Interval Exon позволяют получить экзон (если конец больше начала, то min() - max() и strand=2)
        # Для strand=1 номер координат в Gene Interval Exon прямо позволяют получить экзон
        
        
        # !!!! ГРАНИЦЫ ЭКЗОНОВ ПО ДНК СОВПАДАЮТ ТОЛЬКО С СООТВЕСТВУЮЩЕЙ РЕВИЗВИЕЙ МРНК!!!
        
        rerun=True    
        while (rerun==True):
            try:
                handle = Entrez.efetch(db="nucleotide", id=str(zzz_NC_NAME), rettype="fasta",seq_start=str(from_pos),seq_stop=str(to_pos),strand=strand)
                rerun=False
                dst=np.array(list(SeqIO.FastaIO.SimpleFastaParser(handle)))
                zzz_GENE_SEQ=str(dst[0][1])
            except Exception as e:
                print('error at seq getting for chormosome ' + zzz_NC_NAME)
                time.sleep(20)

        rerun=True        
        while (rerun==True):
            try:
                handle = Entrez.elink(dbfrom="gene",id=str(g_i[count1]),linkname="gene_nuccore_refseqrna")   
                rerun=False
                record =Entrez.read(handle)
            except Exception as e:
                rerun=True
                print('Error at getting mRNA list for ' + g_n[count1])
                time.sleep(20)
        
        
        for R2,count2 in zip(record,range(0,len(g_i))):
            dst=eval(str(R2))
            zzz_mRNA_LIST=[]
            zzz_mRNA_SEQ=[]# список mRNA данного гена
            EXON_B=[]
            CDS_B=[]
            # Через Accession number выходим на идентификаторы
            if (len(dst['LinkSetDb'])>0):
                for i in dst['LinkSetDb'][0]['Link']:
                    rerun=True
                    while (rerun==True):
                        try:
                            handle = Entrez.efetch(db="nucleotide", id=str(i['Id']), rettype="acc")
                            rerun=False
                            xx=handle.readline()
                        except Exception as e:
                            rerun=True
                            print('Error at getting mRNA list for ' + str(i['Id']))
                            time.sleep(20)
                        
            
                    xx=xx.replace('\n','')
                    xx1=xx.split('.');
                    xx1=xx1[0]
                    xx2=xx1.split('_')
                    VOC_LIST.append([g_i[count1],xx1])
                    if ((xx2[0]=='NM') | (xx2[0]=='NR')):
                        zzz_mRNA_LIST.append(xx1)
    
                rerun=True    
                while (rerun==True):
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=zzz_mRNA_LIST, retmode="xml")
                        records_m = Entrez.read(handle)
                        rerun=False
                    except Exception as e:
                        print(e)
                        rerun=True
                        print('Error in download of mRNA seq for ' + g_n_S[count2])
                        time.sleep(20)
                        
                
                for R_m in records_m:
                    FBLOCK=R_m['GBSeq_feature-table']
                    EXON_BLOCK=[]
                    CDS_BLOCK=[]
                    print('Gene: ' + g_n[count1] + ' mRNA: ' + R_m['GBSeq_primary-accession'])
                    for F_m in FBLOCK:
                        if (F_m['GBFeature_key']=='exon'):  
                            st=int(F_m['GBFeature_intervals'][0]['GBInterval_from'])
                            fn=int(F_m['GBFeature_intervals'][0]['GBInterval_to'])
                            pos=[0,0]
                            pos[0]=st
                            pos[1]=fn
                            EXON_BLOCK.append(pos)
                        if (F_m['GBFeature_key']=='CDS'):    
                            st=int(F_m['GBFeature_intervals'][0]['GBInterval_from'])
                            fn=int(F_m['GBFeature_intervals'][0]['GBInterval_to'])
                            pos=[0,0]
                            pos[0]=st
                            pos[1]=fn
                            CDS_BLOCK.append(pos)                            
                    mNAME=R_m["GBSeq_primary-accession"]
                    mDEF=R_m["GBSeq_definition"]
                    zzz_mRNA_SEQ.append(R_m["GBSeq_sequence"])
                    EXON_B.append(EXON_BLOCK)
                    CDS_B.append(CDS_BLOCK)
                    
            else:
                zzz_mRNA_LIST.append('DUMP_RNA')
                zzz_mRNA_SEQ.append(zzz_GENE_SEQ)
                EXON_B.append('DUMP_EXON')
                CDS_B.append('DUMP_CDS')
        
    
            C1=count_list[count1]
            GENE_NAME=g_n[count1]
            GENE_ID=g_i[count1]
            GENE_FUNC=str(f_n[count1])
           
            C2=0
            for Nx,S in zip(zzz_mRNA_LIST,zzz_mRNA_SEQ):
                C2=C2+1   
                TOTAL_C=TOTAL_C+1
                S=S.upper()
                record=SeqRecord(Seq(S, generic_nucleotide),id=Nx,description=str(TOTAL_C) + '_MRNA_' + str(C1) + '_' +  str(C2)+ '_' + GENE_NAME + '_' + str(GENE_ID))
                SeqIO.write(record,'MRNA_FASTA/' + str(C1) + '_' + str(C2) + '.fasta',"fasta")
                TOTAL_M.append(record)
            record=SeqRecord(Seq(zzz_GENE_SEQ, generic_nucleotide),id=str(GENE_NAME),description=str(GENE_ID) + ':' + str(f_n[count1]))
            SeqIO.write(record, 'GENE_FASTA/' + str(C1) + '_' + str(0) + '.fasta', "fasta")
            TOTAL_G.append(record)            
            
            GENE_START=from_pos
            GENE_END=to_pos
            GENE_STRAND=strand
            GENE_STRAND_TEXT=strand_text
            C2=0


            for M,E,C in zip(zzz_mRNA_LIST,EXON_B,CDS_B):
                C2=C2+1
                OUTPUT.append([C1,C2,GENE_ID,GENE_NAME,GENE_FUNC,M,zzz_NC_NAME,GENE_START,GENE_END,GENE_STRAND,GENE_STRAND_TEXT,E,C])


 
    OUTPUT=pd.DataFrame(OUTPUT)           
    OUTPUT=OUTPUT.reset_index(drop=True)
    writer = pd.ExcelWriter(str(CK) + '_SRC_DATA.xlsx')
    OUTPUT.to_excel(writer,index=False)
    writer.save() 
    
    MAIN_BLOCK.append(OUTPUT)
    
MAIN_BLOCK=pd.concat(MAIN_BLOCK)
writer = pd.ExcelWriter('SRC_DATA.xlsx')
MAIN_BLOCK.to_excel(writer,index=False)
writer.save()
VOC=pd.DataFrame(VOC_LIST)
VOC.to_csv('REG_LIST.csv',index=False)