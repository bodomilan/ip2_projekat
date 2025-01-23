from Bio import SeqIO
from typing import List
import numpy as np
import pickle
import pandas as pd
import os
from zipfile import ZipFile
from py7zr import SevenZipFile

# formira fajl koji sadrzi u svakom redu po jednu sekvencu, ako se cita iz CSV fajla vrsi se i preprocesiranje
# trebalo bi da radi i kada se izvrsi ravnanje
def write_sequences_to_file(source_path : str="../data/csv/sars2_mers_sars1.txt", destination_path : str="../data/sequences", is_csv : bool=True):
    sequences = []
    if is_csv:
        X = pd.read_csv(source_path, skipinitialspace = True )
        X.columns = X.columns.str.strip()
    
        sequences = X["KOD_NUC"].tolist()
            
        with open(destination_path, "w") as outfile:
            for seq in sequences:
                outfile.write(seq.strip()+'\n')
    else:
        with open(source_path, "r") as infile:
            for line in infile:
                sequences.append(line)

        with open(destination_path, "w") as outfile:
            for seq in sequences:
                outfile.write(seq)
            

# kao argument dobija tip kompresovanog fajla, fajl sa distancama u obliku trougaone liste listi i kompresuje ga na destination_path
def serialize_and_compress_distance_matrix(source_path: str, destination_path: str, compression_file = ZipFile):
    
    # matrica ce biti u obliku trougla, jer je matrica distanci simetricna pa nema potrebe cuvati je kompletnu
    distances_triangle_list = []
    row_distances = []
    with open(source_path, "r") as infile:
        for line in infile:
            distances = line.split(", ")[:-1]
            for distance in distances:
                row_distances.append(float(distance))
            distances_triangle_list.append(row_distances)
            row_distances = []

    extension = ".zip"
    if(compression_file == SevenZipFile):
        extension = ".7z"
    index = destination_path.find(extension)
    if index == -1:
        index = len(destination_path)
    destination_path = destination_path[:index]
    with open(destination_path, "wb") as outfile:
        pickle.dump(distances_triangle_list, outfile)

    with compression_file(destination_path+extension,'w') as compressor:
        compressor.write(destination_path)

    # obrisi ne kompresovan fajl
    os.remove(destination_path)


# vraca trougaonu listu distanci
def deserialize_and_decompress_distance_matrix(source_path: str) -> List[int]:

    extract_path = ""
    extension_beginning = source_path.rfind(".")
    if(extension_beginning == -1):
        raise Exception("file must contain extension!")
        
    extension = source_path[extension_beginning:]
    if (extension != ".7z" and extension != ".zip"):
        raise Exception("Extension must be .7z or .zip!")

    compression_file = SevenZipFile
    if(extension == ".7z"):
        compression_file = SevenZipFile
        extract_path_index = source_path.rfind("/") + 1
        extract_path = source_path[:extract_path_index]

    elif (extension == ".zip"):
        compression_file = ZipFile
        extract_path_index = source_path.rfind("../")
        if extract_path_index != -1:
            extract_path = source_path[:(extract_path_index+2)]
    
    with compression_file(source_path,'r') as compressor:
        compressor.extractall(extract_path)

    index = source_path.find(extension)
    decompressed_path = source_path[:index]
    with open(decompressed_path, "rb") as infile:
     	distances = pickle.load(infile)

    os.remove(decompressed_path)

    return distances





def five_adic_coding(s: str, cut_len=float('inf')) -> List[int]:
    s.strip() #remove whites from end
    
    seq = []
    num = ''
    i = 0
    for c in s:
        if c == 'C':
            num += '1'
        elif c == 'A':
            num += '2'
        elif c == 'T':
            num += '3'
        else:
            num += '4'
        i += 1
        if i == 3:
            i = 0
            seq.append(int(num))
            num = ''

        if len(seq) == cut_len: #cut_len should be 78
            break
    return seq

def simple_difference_distance(s1: List[int], s2: List[int]) -> int:
    i = len(s1)
    j = 0
    distance = 0

    diff = np.array(s1) - np.array(s2)
    distance = sum([abs(x) for x in diff])
    return distance

def five_two_adic_distance(s1: List[int], s2: List[int]) -> int:
    distance = 0
    duzina1 = len(s1)
    duzina2 = len(s2)

    min_duzina = min(duzina1, duzina2)
    for i in range(min_duzina):
        if int(s1[i]/100) != int(s2[i]/100):
            distance += 1
            
        elif int(s1[i]/10) != int(s2[i]/10):
            distance += 1/5
            
        elif s1[i] % 10 != s2[i] % 10:
            if abs(s1[i] % 10 - s2[i] % 10) == 2:
                distance += 1/2*1/25
            else:
                distance += 1*1/25
                
    return distance + abs(duzina1 - duzina2) / 3.0






# ON NON CSV data
# vraca kodiranu sekvencu i klase proteina
def read_data(filePath: str) -> (List[int], List[str]):
    with open(filePath) as file:
        fasta_sequences_sars1 = SeqIO.parse(file,'fasta')
        coded_sequences = []
        protein_classes = []
        for fasta in fasta_sequences_sars1:
            coded_seq = five_adic_coding(fasta.seq)
            coded_sequences.append(coded_seq)

            protein_class = fasta.description
            # formatiranje proteinske klase da se izvuku samo kljucne informacije o klasi ignorisuci velika i mala slova
            protein_class = ''.join(protein_class.split('|')[1].split('[')[0].strip().split()).lower()
            
            protein_classes.append(protein_class)
    return coded_sequences, protein_classes

def read_dataset_virus(file_path_sars1="../data/sars1.fasta", file_path_mers="../data/mers.fasta",  file_path_sars2="../data/sars2.fasta"):
    """
        Vraca niz, vraca y labele (tip virusa) i y_protein labele (tip proteina)
    """
    X_sars1, y_protein_sars1 = read_data(file_path_sars1)
    y_sars1 = ['sars1' for i in range(len(X_sars1))]
    X_mers, y_protein_mers = read_data(file_path_mers)
    X_sars2, y_protein_sars2 = read_data(file_path_sars2)
    
    X=X_sars1
    y_protein = y_protein_sars1
    y=y_sars1
    
    for i in range(len(X_mers)):
        X.append(X_mers[i])
        y_protein.append(y_protein_mers[i])
        y.append('mers')
    
    for i in range(len(X_sars2)):
        X.append(X_sars2[i])
        y_protein.append(y_protein_sars2[i])
        y.append('sars2')

    
    return X,y,y_protein

def train_test_virus():
    y_train=pd.read_csv("../data/virus/y_train.csv",index_col=0)
    y_test=pd.read_csv("../data/virus/y_test.csv",index_col=0)
    X_train,X_test=[],[]
    with open("../data/virus/X_train.pkl", "rb") as file:
        X_train=pickle.load(file)
    with open("../data/virus/X_test.pkl", "rb") as file:
        X_test=pickle.load(file)
    return X_train, X_test, y_train, y_test

def train_test_protein():
    y_train=pd.read_csv("../data/protein/y_train.csv",index_col=0)
    y_test=pd.read_csv("../data/protein/y_test.csv",index_col=0)
    X_train,X_test=[],[]
    with open("../data/protein/X_train.pkl", "rb") as file:
        X_train=pickle.load(file)
    with open("../data/protein/X_test.pkl", "rb") as file:
        X_test=pickle.load(file)
    return X_train, X_test, y_train, y_test

def train_test_sars2():
    y_train=pd.read_csv("../data/sars2_who/y_train.csv",index_col=0)
    y_test=pd.read_csv("../data/sars2_who/y_test.csv",index_col=0)
    X_train,X_test=[],[]
    with open("../data/sars2_who/X_train.pkl", "rb") as file:
        X_train=pickle.load(file)
    with open("../data/sars2_who/X_test.pkl", "rb") as file:
        X_test=pickle.load(file)
    return X_train, X_test, y_train, y_test