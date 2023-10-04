
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#Fungsi menghitung basa nitrogen 
def count_each_nitrogen_base(genome_sequence):
    Adenin = 0
    Guanin = 0
    Cytosin = 0
    Timin = 0 
    Error = 0
    for base in genome_sequence:
        if base == "A":
            Adenin+=1
        elif base == "G":
            Guanin+=1
        elif base == "C":
            Cytosin+=1
        elif base == "T":
            Timin+=1
        else:
            Error+= 1
    return Adenin, Guanin, Cytosin, Timin, Error


def count_dimer(window, genome_sequence):
    window_sequence = []
    current_string = ""
    for base in genome_sequence:
        current_string += base
        if len(current_string) == window:
            window_sequence.append(current_string)
            current_string = ""
    window_sequence.append(current_string)
    gc= np.zeros(len(window_sequence))
    ta= np.zeros(len(window_sequence))
    for piece,sequence in enumerate(window_sequence):
        for base in sequence:
            if base=="G" or base =="C":
                gc[piece] += 1
            elif base=="T" or base =="A":
                ta[piece] += 1

    gc_ta_frequency_in_windo_sequence = {}
    for i in range(len(window_sequence)):
        window_info = {
            "Sequence": window_sequence[i],
            "Number GC": gc[i]/window,
            "Number TA": ta[i]/window
        }
        gc_ta_frequency_in_windo_sequence[f"Window_{i}"] = window_info
    return gc_ta_frequency_in_windo_sequence



# Example usage:
file_path = "C:\Database\Semester 5\Komputasi Genomik\H_influenzae.fasta"  # Replace with the path to your FASTA file
#genome_data = read_fasta(file_path)

genome_sequence = ""
with open(file_path, 'r') as file:
    for line in file:
        if not line.startswith(">"):  # Skip header lines starting with '>'
            genome_sequence += line.strip()
        else:
            desc = line.strip()
# print("Genome Data:")
print(desc)
#print(genome_sequence)

Adenin, Guanin, Cytosin, Timin, Error = count_each_nitrogen_base(genome_sequence)

print(len(genome_sequence))
print(f"Jumlah Adenin: {Adenin} basa")
print(f"Jumlah Guanin: {Guanin} basa")
print(f"Jumlah Cytosin: {Cytosin} basa")
print(f"Jumlah Timin: {Timin} basa")
print(f"Error rate: {Error/(len(genome_sequence))} %")

#Moving Window Length
window = int(input("Panjang wwindow yang diinginkan:"))
gc_ta_frequency_in_windo_sequence= count_dimer(window, genome_sequence)

# Menyiapkan data dari kamus gc_ta_frequency_in_windo_sequence
data = []
for window_key, window_info in gc_ta_frequency_in_windo_sequence.items():
    data.append([window_key, window_info["Sequence"], window_info["Number GC"], window_info["Number TA"]])

# Membuat DataFrame
data_dimer = pd.DataFrame(data, columns=['Window', 'Sequence', 'Number GC', 'Number TA'])

# Menampilkan DataFrame
print(data_dimer)

# Plotting GC and TA frequencies
x = range(len(gc_ta_frequency_in_windo_sequence))
gc_values = [item["Number GC"] for item in gc_ta_frequency_in_windo_sequence.values()]
ta_values = [item["Number TA"] for item in gc_ta_frequency_in_windo_sequence.values()]

plt.figure(figsize=(10, 5))
plt.plot(x, gc_values, label="GC Frequency")
plt.plot(x, ta_values, label="TA Frequency")
plt.xlabel("Window Index")
plt.ylabel("Frequency")
plt.title("GC and TA Frequency in Windows")
plt.legend()
plt.grid()
plt.show()





    




    



