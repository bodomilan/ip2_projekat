from Bio import SeqIO
import numpy as np

def five_adic_coding(s):
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

        if len(seq) == 3: # should be 40
            break
    return seq

def read_data(filePath):
    with open(filePath) as file:
        fasta_sequences_sars1 = SeqIO.parse(file,'fasta')
        data = []
        for fasta in fasta_sequences_sars1:
            seq = five_adic_coding(fasta.seq)
            data.append(seq)
    return data

def simple_difference_distance(s1, s2):
    i = len(s1)
    j = 0
    distance = 0

    diff = np.array(s1) - np.array(s2)
    distance = sum([abs(x) for x in diff])
    return distance

def five_two_adic_distance(s1, s2):
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

def read_dataset_virus(file_path_sars1="../data/sars1.fasta", file_path_mers="../data/mers.fasta",  file_path_sars2="./data/sars2.fasta"):
    """
        Vraca niz 
    """
    X_sars1 = read_data("../data/sars1.fasta")
    y_sars1 = ['sars1' for i in range(len(X_sars1))]
    X_mers = read_data("../data/mers.fasta")
    #X_sars2=fdp.read_data("./data/sars2.fasta")
    X=X_sars1
    y=y_sars1
    for x in X_mers:
        X.append(x)
        y.append('mers')
    
#     for x in X_sars2:
#         X.append(x)
#         y.append('sars2')

    
    return X,y

