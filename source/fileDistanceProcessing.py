from Bio import SeqIO
from typing import List
import numpy as np

def five_adic_coding(s: str) -> List[int]:
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

        if len(seq) == 39: # should be 39
            break
    return seq

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

def simple_difference_distance(s1: List[int], s2: List[int]) -> int:
    i = len(s1)
    j = 0
    distance = 0

    diff = np.array(s1) - np.array(s2)
    distance = sum([abs(x) for x in diff])
    return distance

def five_two_adic_distance(s1: List[int], s2: List[int]) -> int:
    distance = 0
    for i in range(len(s1)):
        if int(s1[i]/100) != int(s2[i]/100):
            distance += 1
            
        elif int(s1[i]/10) != int(s2[i]/10):
            distance += 1/5
            
        elif s1[i] % 10 != s2[i] % 10:
            if abs(s1[i] % 10 - s2[i] % 10) == 2:
                distance += 1/2*1/25
            else:
                distance += 1*1/25
    return distance

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

