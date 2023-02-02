#!/usr/bin/env python
# coding: utf-8

# In[69]:


from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW
import requests
import re
from Bio import SeqIO
from Bio import Entrez
from Bio import Medline
from tqdm import tqdm
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from bs4 import BeautifulSoup
from Bio.SeqFeature import CompoundLocation
import mysql.connector as SQLC
import mysql.connector


# In[70]:


DataBase = SQLC.connect(
    host ="geo.di.uminho.pt",
    user ="bioinformatica",
    password ="20221207",
    database ="AP_db_KRG"
)
DataBase.autocommit = True 
Cursor = DataBase.cursor()

num= input("Escolha o tamanho máximo das seq: ")
Seq_and_length=[]
just_seq=[]
try:
    sql_Seq_length= f"select Gene.ID_genebank, Gene.length from Gene where Gene.length <{num}"
    Cursor.execute(sql_Seq_length)
    for row in Cursor:
        Seq_and_length.append(str(row)) 
    
    sql_Seq= f"select Gene.ID_genebank from Gene where Gene.length <{num}"
    Cursor.execute(sql_Seq)
    for row in Cursor:
        just_seq.append(str(row))
        
except mysql.connector.Error as e:
    print("Erro na escrita na base de dados: {}".format(e) )    
    
print (Seq_and_length)

if just_seq == []:
        print()
        print('nenhuma seq inferior a esse tamanho, pesquise de novo.')


# In[71]:


list_seq=[]
choose= input('escolha a primeira seq com que quer trabalhar: ')
inter= int(choose)
extr= ', '.join(just_seq)
h= extr.replace("(",'')
hh= h.replace(",)",'')
ID_al= hh.split(', ')
select= ID_al[inter]
list_seq.append(select)

# choose2= input('escolha a segunda seq com que quer trabalhar: ')
# inter2= int(choose2)
# extr2= ', '.join(just_seq)
# h2= extr2.replace("(",'')
# hh2= h2.replace(",)",'')
# ID_al2= hh2.split(', ')
# select2= ID_al2[inter2]
# list_seq.append(select2)

limpar_aspas_all= ', '.join(ID_al)
retirar_all= limpar_aspas_all.replace("'","")
otput_final_all= retirar_all.split(', ')

limpar_aspas= ', '.join(list_seq)
retirar= limpar_aspas.replace("'","")
otput_final= retirar.split(', ')


#cria lista com todos as seq com o tamanho máximo definido:
todas_as_seq=[]
database = 'nucleotide'
email= 'rodrigoce9@gmail.com'
idlist= otput_final_all
handle = Entrez.efetch(db=database, id=idlist, rettype="gb") 
records = list(SeqIO.parse(handle,"gb"))
handle.close()
for info in records:
    seqs= (info.seq)
    todas_as_seq.append(seqs)
    
choose2= input('tipo de molécula com que está a trabalhar: ')    
#cria lista só com as duas seq escolhidas para o alinhamento 
seq_for_al=[]
database = 'nucleotide'
email= 'rodrigoce9@gmail.com'
idlist1= otput_final[0]
# idlist2= otput_final[1]
handle1 = Entrez.efetch(db=database, id=idlist1, rettype="gb") 
records1 = SeqIO.read(handle1,"gb")
handle1.close()
s1_=records1.seq
# handle2 = Entrez.efetch(db=database, id=idlist2, rettype="gb") 
# records2 = SeqIO.read(handle2,"gb")
# handle2.close()
# s2_=records2.seq
print(todas_as_seq)
print(s1_)
# print(s2_)


# __Goal__ Implementation of sequence manipulation algorithms: transcription, translation, ORFs, extraction of all potential proteins

# In[72]:


#function that allows us to associate a codon with an amino acid- Human
def translate_codon(cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", 
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod] 
    else: return "Don´t exist"
translate_codon("xbt")


# In[73]:


#talvez colocar mais possibilidade porque eu penso que a tradução muda para cada organismo??


# In[96]:


#a class that allows to accomplish the objective
class sequence_info:
    #"building block", where you always have to have the sequence you are analyzing and the type of molecule: dna, rna, protein
    def __init__(self,sequence:s1_, sequence_type:choose2):
        if sequence_type not in ["DNA","RNA","PROTEIN"]:
            raise ValueError("Invalid sequence type")
        self.sequence=sequence.upper() #uppercase sequence as a universal rule
        self.sequence_type=sequence_type
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self,n):
        return self.sequence[n]
    def __getslice__(self,i,j):
        return self.sequence[i:j]
    def __str__(self):
        return self.sequence
    def sequence_type(self):
        return self.sequence_type
    def letters(self):
        #each letter possible to have for dna, rna , protein for later to be able to validate the information about the type of molecule corresponding to the sequence
        if (self.sequence_type=="DNA"): return "ACGT"
        elif(self.sequence_type=="RNA"):return "ACGU"
        elif(self.sequence_type=="PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
        else: raise ValueError("ERROR")
    def validate(self): #validates information about the type of molecule associated with the sequence
        letras=self.letters()
        res=True
        i=0
        while i<len(self.sequence) and res:
            if self.sequence[i] not in letras: res = False
            else: i +=1
        if res==False:
            raise ValueError("Put the right name of a molecule")
        else: return res 
    def contagem(self):
        contagem={}
        for x in self.sequence:
            contagem[x]=contagem.get(x,0)+1
        return sorted(contagem.items(), key=lambda x: x[1], reverse=True)      
    def reverse_complement(self):#falta comprovar isto
        if (self.validate() == False) or (self.sequence_type in ["RNA","PROTEIN"]) :return "only do reverse complement with DNA"
        t=self.sequence.replace("T","a")
        g=t.replace("G","c")
        c=g.replace("C","g")
        a=c.replace("A","t")
        reverse_comp=a[::-1]
        return reverse_comp.upper()    
    def transcription(self):
        if (self.validate() == False) or (self.sequence_type != "DNA"):return "only do transcription with DNA"
        else: return  sequence_info(self.sequence.replace("T","U"),"RNA")
    def translate(self):
        if (self.validate() == False) or (self.sequence_type != "DNA"): return "only want to do translation with DNA"
        sequence_aa=""
        for x in range(0, len(self.sequence)-2,3):#3 by 3
            codon=self.sequence[x:x+3]
            sequence_aa += translate_codon(codon)
        return sequence_aa
    def get_orfs_DNA(self): #6 possible orfs
        if (self.validate() == False) or (self.sequence_type != "DNA"): return "only want DNA orf"
        inverse_sequence=self.sequence[::-1]
        return ",".join([self.sequence[p:]for p in range(3)] + [inverse_sequence[p:] for p in range(3)]) 
    def get_orfs_Proteins(self): # translation of the 6 possible dna orfs
        if (self.validate() == False) or (self.sequence_type != "DNA"): return "just work from dna orf"
        inverse_sequence=self.sequence[::-1]
        orfs_protein=[]
        for p in range(3):
            orfs_protein.append(sequence_info(self.sequence[p:],self.sequence_type).translate())
            orfs_protein.append(sequence_info(inverse_sequence[p:],self.sequence_type).translate())
        return ",".join(orfs_protein)
    def protsInAA(self):#gives all possible proteins of an amino acid sequence( from M to _)
        seq_aa = self.translate()
        prots_act = []
        prots = []
        for aa in seq_aa:
            if aa == '_':
                if prots_act:
                    for p in prots_act:
                        prots.append(p)
                    prots_act = []
            else:
                if aa == 'M':
                    prots_act.append('')
                for i in range(len(prots_act)):
                    prots_act[i] += aa
        return ",".join(prots) 
 #if name=main é para conseguir importar direitinho os nodulos ou seja quanto estou neste ficheiro o name =main e ele faz a classe se for para outro ficheiro e importar a classe ele faz na mesma a classe sem repetir tufo o que aqui tenho
# https://www.youtube.com/watch?v=3dHyS1W4TIE
if __name__ == "__main__":
    s1_="AA??"
    s1 = sequence_info(s1_,choose2)
    s2 = sequence_info("MKVVLSVQERSVVSLL_", "PROTEIN")
    #s3 = sequence_info("ATTTTBTT","PROTEIN")
    #print("sequence 1 is {}? {}".format(s1.sequence_type,s1.validate()))
    s3 = s1.validate()
    print(s3)
    #s4 = s1.reverse_complement()
    #print(s4)
#     s5= s1.translate()
#     print(s5)
#     s6=s1.get_orfs_DNA()
#     print(s6)
#     s8=s1.get_orfs_Proteins()
#     print(s8)
#     s9=s1.protsInAA()
#     print(s9)
#     s10=s1.contagem()
#     print(s10)
    #s11=s1.__len__()
    #print(s11)
    s12=s1.letters()
    print(s12)
    

    
    


# In[ ]:




