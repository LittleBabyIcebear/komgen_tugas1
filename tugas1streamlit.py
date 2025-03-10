import streamlit as st
import numpy as np
import pandas as pd
import random
import plotly.express as px
import plotly.graph_objects as go
from Basanitrogen import Codon_DNA
import os
import shutil

def generate_random_sequence(sequence_type, nitrogen_base, num_sequences):
    random_sequence = []
    for i in range(0, num_sequences):
        number = random.randint(0, len(nitrogen_base) - 1)
        if sequence_type == "RNA" and number == 3:
            number = 4
        elif sequence_type == "DNA" and number == 4:
            number = 3
        random_sequence.append(nitrogen_base[number])
    return random_sequence
    
def split_sequence_into_codons(random_sequence):
    split_sequence = []
    current_string = ""
    for base in random_sequence:
        current_string += base
        if len(current_string) == 3:
            split_sequence.append(current_string)
            current_string = ""
    return split_sequence

def map_codons_to_amino_acids(split_sequence, Codon_DNA):
    for base in split_sequence:
        for amino_acid, data in Codon_DNA.items():
            if base in data["Codon"]:
                st.write(f"Codon {base}: {amino_acid} ({data['Single_Letter']})")

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
    a = np.zeros(len(window_sequence))
    g = np.zeros(len(window_sequence))
    c = np.zeros(len(window_sequence))
    t = np.zeros(len(window_sequence))
    for piece,sequence in enumerate(window_sequence):
        for base in sequence:
            if base=="G":
                g[piece]+=1 
                gc[piece] += 1
            elif base =="C":
                c[piece]+= 1
                gc[piece] += 1
            elif base=="T":
                t[piece] += 1
                ta[piece] += 1                
            elif base =="A":
                a[piece] += 1
                ta[piece] += 1

    gc_ta_frequency_in_windo_sequence = {}
    for i in range(len(window_sequence)):
        window_info = {
            "Sequence": window_sequence[i],
            "Number GC": gc[i],
            "Number TA": ta[i],
            "Freq GC": gc[i]/(window),
            "Freq TA": ta[i]/(window),
            "Freq A": a[i]/window,
            "Freq G": g[i]/window,
            "Freq T": c[i]/window,
            "Freq C": t[i]/window,
        }
        gc_ta_frequency_in_windo_sequence[f"Window_{i}"] = window_info
    return gc_ta_frequency_in_windo_sequence

def change_point_analysis(gc_values, ta_values, gc_ta_frequency_in_windo_sequence):
    # Initial Condition
    change_point = [0] # Inisialisasi change_point dengan nilai awal 0
    change_point.append(0)
    if gc_values[0] > ta_values[0]:
        top = gc_values
        bottom = ta_values

    elif gc_values[0] < ta_values[0]:
        top = ta_values
        bottom = gc_values

    for i in range(0, len(gc_ta_frequency_in_windo_sequence)):
        if top[i] > bottom[i] and change_point[i-1] != 1:
            change_point.append(1)
        elif top[i] < bottom[i] and change_point[i-1] != 0:
            change_point.append(0)
        else:
            change_point.append(change_point[i-1])
    return change_point

def orf_finder(genome_sequence, Codon_DNA, orf_code, index_orf, min_length, plus_or_min):
    def translate_sequence_into_codons(sequence,Codon_DNA):
        split_sequence = []
        sequence_single_letter=[]
        current_string = ""
        for base in sequence:
            current_string += base
            if len(current_string) == 3:
                split_sequence.append(current_string)
                current_string = ""
        for base in split_sequence:
            for amino_acid, data in Codon_DNA.items():
                if base in data["Codon"]:
                    sequence_single_letter.append(data['Single_Letter'])
        return sequence_single_letter

    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    final_dictionary = {}
    index_orf_in_each_code = 1
    
    if plus_or_min == "+": 
        operator_frame = "+"
    elif plus_or_min == "-":
        operator_frame = "-"

    if orf_code == 1:
        i = 0
        type_frame = 1
    elif orf_code == 2:
        i = 1
        type_frame = 2
    elif orf_code == 3:
        i = 2
        type_frame = 3
    
    index_orf += 1
    
    while i < len(genome_sequence):
        if i + 2 < len(genome_sequence) and genome_sequence[i] == 'A' and genome_sequence[i+1] == 'T' and genome_sequence[i+2] == 'G':
            j = i + 3
            while j < len(genome_sequence) - 2:
                codon = genome_sequence[j] + genome_sequence[j+1] + genome_sequence[j+2]
                if codon in stop_codons:
                    orf_sequence = genome_sequence[i:j+3]
                    if len(orf_sequence) >= min_length:
                        if operator_frame == "+":
                            orf_info = {
                                "Index_ORF" : index_orf,
                                "Frame ORF": f" Frame [{operator_frame}{type_frame}] Number in Each Frame: {index_orf_in_each_code}",
                                "Index Start Codon": i+1,
                                "Index Stop Codon": j+3,
                                "ORF Lenght" :(j+2) - (i+1) +2,
                                "Sequence ORF": orf_sequence,
                                "Sequence Amino Acid": "".join(translate_sequence_into_codons(orf_sequence, Codon_DNA))
                            }
                        elif operator_frame == "-":
                            orf_info = {
                                "Index_ORF" : index_orf,
                                "Frame ORF": f" Frame [{operator_frame}{type_frame}] Number in Each Frame: {index_orf_in_each_code}",
                                "Index Start Codon": len(genome_sequence) - i,
                                "Index Stop Codon": len(genome_sequence) - j-2,
                                "ORF Lenght" :(j+2) - (i+1) +2,
                                "Sequence ORF": orf_sequence,
                                "Sequence Amino Acid": "".join(translate_sequence_into_codons(orf_sequence, Codon_DNA))
                            }
                        final_dictionary[f"Index ORF_{index_orf_in_each_code}"] = orf_info
                        index_orf += 1
                        index_orf_in_each_code += 1
                    i = j + 3
                    break
                j += 3
            else:
                i += 3
        else:
            i += 3
    
    return final_dictionary

def main():
    purines = ["A", "G"]
    pyrimidines = ["C", "T", "U"]
    nitrogen_base = purines + pyrimidines
    st.sidebar.title("Choose It! 🔎")

    selected_option = st.sidebar.radio("Select an option", ["Show Amino Acid Table", "Analyze Sequence to Amino Acid", "Analyze Fasta File", "Find Open Reading Frame"])

    if selected_option == "Show Amino Acid Table":
        st.title("🧬Random DNA/RNA Sequence Generator🧬")
        st.subheader("Nitrogen Base Classification")
        st.write(f"Purines: {purines}")
        st.write(f"Pyrimidines: {pyrimidines}")
        st.write("DNA: Adenin, Gunanin, Cytosin, Timin")
        st.write("RNA: Adenin, Gunanin, Cytosin, Uracil")
        st.subheader("Related Codon to Amino Acid Table")
        df = pd.DataFrame([(amino_acid, ', '.join(data['Codon']), data['Single_Letter']) for amino_acid, data in Codon_DNA.items()],
                            columns=['Amino Acid', 'Codon', 'Single Letter'])

        st.table(df)

    elif selected_option == "Analyze Sequence to Amino Acid":
        st.title("🧬Random DNA/RNA Sequence Generator🧬")
        choose = st.radio("Choose RNA or DNA sequence?", ("RNA", "DNA"))
        num_sequences = st.number_input("Enter the number of sequences (must be a multiple of 3):", min_value=3, step=3)
        random_sequence = generate_random_sequence(choose, nitrogen_base, num_sequences)

        #result_container = st.empty()  # Wadah untuk hasil saat sudah diubah menjadi 

        if st.button("Run"):
            st.subheader("Generated Sequence")
            st.write("".join(random_sequence))

            if choose == "RNA":
                random_sequence = [base if base != "U" else "T" for base in random_sequence]
            
            split_sequence = split_sequence_into_codons(random_sequence)

            st.subheader("Split Sequence into Codons")
            st.write(split_sequence)

            st.subheader("Codon to Amino Acid Mapping")
            map_codons_to_amino_acids(split_sequence, Codon_DNA)
    
    elif selected_option == "Analyze Fasta File":
        st.title("🧬Adventure to the Fasta File🧬")
        st.subheader("Uploaded FASTA File")
        #Elemen upload file
        uploaded_file = st.file_uploader("Drag FASTA file here!", type=["fasta", "fa"])
        window = st.number_input("Given window length", min_value=1)
        if st.button("Run"):
            genome_data_with_header = uploaded_file.read().decode()

            #Deskripsi untuk line pertama  
            st.subheader("Description of the Genom")
            lines = genome_data_with_header.split("\n")
            if lines:
                desc = lines[0]
                st.write(desc)
            
            genome_data_no_header = "".join(genome_data_with_header.split("\n")[1:])  # Remove the first line (header)
            Adenin, Guanin, Cytosin, Timin, Error = count_each_nitrogen_base(genome_data_no_header)
            
            st.subheader("Genome Metrics")
            # Create a pie chart for Adenin, Guanin, Timin, and Cytosin using Plotly
            st.write(f"Length of the data: {len(genome_data_no_header)} base pair (bp)")
            data = {
                 "Parameter": ["Adenin", "Guanin", "Cytosin", "Timin", "Error rate"],
                 "Value (base)": [Adenin/ len(genome_data_no_header), Guanin/ len(genome_data_no_header), Cytosin/ len(genome_data_no_header), Timin/ len(genome_data_no_header), Error / len(genome_data_no_header)]
            }
            df = pd.DataFrame(data)
            st.table(df)
            
            fig = px.pie(
                names=['Adenin', 'Guanin', 'Timin', 'Cytosin'],
                values=[Adenin, Guanin, Timin, Cytosin],
                title='Genome Metrics: Adenin, Guanin, Timin, and Cytosin'
            )

            # Display the pie chart using Streamlit
            st.plotly_chart(fig)

            # Tampilan judul
            st.title("Analyze GC TA Dimers in Window Genom")

            gc_ta_frequency_in_windo_sequence = count_dimer(window, genome_data_no_header)
            
            # # Menyiapkan data dari kamus gc_ta_frequency_in_windo_sequence
            data = []
            
            # Membuat DataFrame
            data_dimer = pd.DataFrame(data, columns=['Window', 'Sequence', 'Number GC', 'Number TA'])

            # Menampilkan DataFrame
            
            x = range(len(gc_ta_frequency_in_windo_sequence))
            gc_values = [item["Freq GC"] for item in gc_ta_frequency_in_windo_sequence.values()]
            ta_values = [item["Freq TA"] for item in gc_ta_frequency_in_windo_sequence.values()]
            a_values = [item["Freq A"] for item in gc_ta_frequency_in_windo_sequence.values()]
            g_values = [item["Freq G"] for item in gc_ta_frequency_in_windo_sequence.values()]
            c_values = [item["Freq C"] for item in gc_ta_frequency_in_windo_sequence.values()]
            t_values = [item["Freq T"] for item in gc_ta_frequency_in_windo_sequence.values()]

            # Membuat DataFrame dari data
            data_monomer = pd.DataFrame({'x': x,'G Frequency': g_values, 'C Frequency': c_values, 'T Frequency': t_values, 'A Frequency': a_values} )

            # Membuat plot Plotly
            fig_monomer = px.line(data_monomer, x='x', y=['G Frequency', 'C Frequency', 'T Frequency', 'A Frequency'], title="Frequency Graphic Monomer A, G, C and T")

            # Menambahkan label pada sumbu x dan sumbu y
            fig_monomer.update_xaxes(title_text="Window")
            fig_monomer.update_yaxes(title_text="Frequency")
            fig_monomer.update_layout(xaxis=dict(fixedrange=False), yaxis=dict(fixedrange=False))  # Mengaktifkan fitur zoomable
            st.plotly_chart(fig_monomer, use_container_width=True) # Menampilkan plot Plotly di Streamlit

            # Membuat DataFrame dari data
            data = pd.DataFrame({'x': x, 'GC Frequency': gc_values, 'TA Frequency': ta_values})

            # Membuat plot Plotly
            fig = px.line(data, x='x', y=['GC Frequency', 'TA Frequency'], title="Frequency Graphic Dimer GC and TA")

            # Menambahkan label pada sumbu x dan sumbu y
            fig.update_xaxes(title_text="Window")
            fig.update_yaxes(title_text="Frequency")

            # Mengaktifkan fitur zoomable
            fig.update_layout(xaxis=dict(fixedrange=False), yaxis=dict(fixedrange=False))

            # Menampilkan plot Plotly di Streamlit
            st.plotly_chart(fig, use_container_width=True)

            #Change Point Analysis
            st.subheader("Change Point When Freq GC and TA Change")
            change_point = change_point_analysis(gc_values, ta_values, gc_values)
            fig_change_point = go.Figure()
            fig_change_point.add_trace(go.Scatter(x=list(x), y=gc_values, mode='lines', name="GC Frequency"))
            fig_change_point.add_trace(go.Scatter(x=list(x), y=ta_values, mode='lines', name="TA Frequency"))
            fig_change_point.add_trace(go.Scatter(x=list(x), y=change_point, mode='lines+markers', name="Change Point"))
            fig_change_point.update_layout(
                xaxis_title="Window Index",
                yaxis_title="Frequency",
                title="GC and TA Frequency in Windows",
                showlegend=True,
                xaxis=dict(showgrid=True),
                yaxis=dict(showgrid=True)
            )

            st.plotly_chart(fig_change_point, use_container_width=True)

    elif selected_option == "Find Open Reading Frame":
        st.title("🧬Find Gene in the DNA Sequence🧬")
        st.subheader("Uploaded FASTA File")
        #Elemen upload file
        uploaded_file = st.file_uploader("Drag FASTA file here!", type=["fasta", "fa"])
        min_length = st.number_input("Given minimal lenght", min_value=1)
        st.write("k_value is a variabel that define the nitrogen base minimum lenght of the ORF")
        if st.button("Run"):
            genome_data_with_header = uploaded_file.read().decode()

            #Deskripsi untuk line pertama  
            st.subheader("Description of the Genom")
            lines = genome_data_with_header.split("\n")
            if lines:
                desc = lines[0]
                st.write(desc)
            
            genome_sequence = "".join(genome_data_with_header.split("\n")[1:])  
            # Membuat urutan baru untuk hasil pembalikan
            inverse_genome_sequence = ""

            # Melakukan pembalikan urutan dan penggantian nukleotida
            for nucleotide in reversed(genome_sequence):
                if nucleotide == "A":
                    inverse_genome_sequence += "T"
                elif nucleotide == "T":
                    inverse_genome_sequence += "A"
                elif nucleotide == "C":
                    inverse_genome_sequence += "G"
                elif nucleotide == "G":
                    inverse_genome_sequence += "C"

            orf_positif_1 = orf_finder(genome_sequence, Codon_DNA, 1, 0, min_length, "+")
            orf_positif_2 = orf_finder(genome_sequence, Codon_DNA, 2, len(orf_positif_1), min_length, "+")
            orf_positif_3 = orf_finder(genome_sequence, Codon_DNA, 3, len(orf_positif_2) + len(orf_positif_1), min_length, "+")
            orf_negatif_4 = orf_finder(inverse_genome_sequence, Codon_DNA, 1, len(orf_positif_2) + len(orf_positif_1)+ len(orf_positif_3), min_length, "-")
            orf_negatif_5 = orf_finder(inverse_genome_sequence, Codon_DNA, 2, len(orf_negatif_4)+len(orf_positif_2) + len(orf_positif_1)+ len(orf_positif_3), min_length, "-")
            orf_negatif_6 = orf_finder(inverse_genome_sequence, Codon_DNA, 3, len(orf_negatif_5) + len(orf_negatif_4)+len(orf_positif_2) + len(orf_positif_1)+ len(orf_positif_3), min_length, "-")
            combined_dicts = [orf_positif_1, orf_positif_2, orf_positif_3, orf_negatif_4, orf_negatif_5, orf_negatif_6]

            st.subheader("Number ORF Found")

            st.write(f"(ORF +1):   {len(orf_positif_1)}")
            st.write(f"(ORF +2):   {len(orf_positif_2)}")
            st.write(f"(ORF +3):   {len(orf_positif_3)}")
            st.write(f"(ORF -1):   {len(orf_negatif_4)}")
            st.write(f"(ORF -2):   {len(orf_negatif_5)}")
            st.write(f"(ORF -3):   {len(orf_negatif_6)}")
            st.write(f"Total ORF Found:   {len(orf_positif_1)+len(orf_positif_2)+len(orf_positif_3)+len(orf_negatif_4)+len(orf_negatif_5)+len(orf_negatif_6)}")

            st.subheader("ORF List:")
            for index, dictionary in enumerate(combined_dicts, start=1):
                st.write("--------")
                for key, value in dictionary.items():
                    if isinstance(value, dict):
                        st.text("  " + key + ":")
                        for inner_key, inner_value in value.items():
                            if inner_key in ["Frame ORF", "Index Start Codon", "Index Stop Codon", "ORF Lenght"]:
                                st.markdown(f"    **{inner_key}:** {inner_value}")
                            else:
                                st.markdown(f"    **{inner_key}:** {inner_value}")
                    else:
                        st.text("  " + key + ": " + value)
            

# Fungsi untuk mengunduh data sebagai CSV
        def download_combined_dict(data, file_name):
            df = pd.DataFrame(data)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Data as CSV",
                data=csv,
                key=file_name,
                file_name=file_name,
                mime='text/csv'
            )

        download_combined_dict(combined_dicts, "combined_dict.csv")

    
if __name__ == "__main__":
    main()
