import os 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.PDBParser import PDBParser
from ramachandraw.parser import get_phi_psi
from ramachandraw.utils import fetch_pdb, plot

#os things
pdb_entry = input("Wprowadź PDB ID białka które chcesz zbadać: ")
os.mkdir(pdb_entry) #tworzy folder o nazwie pdb_entry
os.chdir(pdb_entry) #zmienia PWD na folder pdb_entry, wszystkie dane będą zapisywane w tym folderze

# pierwsza funkcja, pobiera PDB
def Download(pdb_entry):
    url = f"https://files.rcsb.org/download/{pdb_entry}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_entry}.pdb", "wb") as file:
            file.write(response.content)
            print(f"Plik został pobrany pomyślnie w folderze {pdb_entry}!")
    else:
        print("Nie udało się pobrac pliku!")
Download(pdb_entry)

# druga funkcja zapisuje sekwencję białka do pliku seq_pdb_entry.txt
def GetSequence(pdb_entry):
    sequence = f'{pdb_entry}.pdb'
    for sequence in SeqIO.parse(sequence, "pdb-seqres"):
        file1 = open(f'seq_{pdb_entry}.txt', 'w')
        file1.write(f'Sekwencja białka {pdb_entry}:\n {sequence.seq}')
        file1.close()
        print(f'Sekwencja została pomyślnie zapisana do pliku seq_{pdb_entry}.txt!')

GetSequence(pdb_entry)

#trzecia funkcja liczy częstotliwość pojawiania się aminokwasu w sekwencji i zapisuje wykres do pliku Plot1.png
def AAcount(pdb_entry):
    sequence = f'{pdb_entry}.pdb'
    for sequence in SeqIO.parse(sequence, "pdb-seqres"):
        sequence = str(sequence.seq)
        analysed_sequence = ProteinAnalysis(sequence)
        aacount = analysed_sequence.count_amino_acids()
        labels = list(aacount.keys())
        sizes = list(aacount.values())
        fig, ax = plt.subplots(figsize=(10,10))
        ax.pie(sizes,labels=labels,autopct='%1.1f%%', pctdistance=1.2, labeldistance=0.6, startangle=90)
        plt.title(f'Udział procentowy aminokwasów w sekwencji białka {pdb_entry}')
        plt.savefig(f'Sequence_plot_{pdb_entry}.png')
        print(f'Procentowa zawartość aminokwasów została pomyślnie zapisana na wykresie pliku Sequence_plot_{pdb_entry}.png!')

AAcount(pdb_entry)

#czwarta funkcja rysuje wykres Ramachandrana na podstawie danych w PDB
def RamachandranPlot(pdb_entry):
    ram_plot = f"{pdb_entry}.pdb"
    plot(ram_plot, save=True, show=False, filename=f'Ramachandran_plot_{pdb_entry}.png')
    print(f'Wyres Ramachandrana został pomyślnie zapisany do pliku Ramachandran_plot_{pdb_entry}.png!')

RamachandranPlot(pdb_entry)