import streamlit as st
import pandas as pd
import random
import plotly.express as px
from Basanitrogen import Codon_DNA

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


def main():
    purines = ["A", "G"]
    pyrimidines = ["C", "T", "U"]
    nitrogen_base = purines + pyrimidines
    st.title("🧬Random DNA/RNA Sequence Generator🧬")
    st.sidebar.title("Choose It! 🔎")

    selected_option = st.sidebar.radio("Select an option", ["Show Amino Acid Table", "Analyze Sequence to Amino Acid", "Analyze Fasta File"])

    if selected_option == "Show Amino Acid Table":
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
        st.subheader("Uploaded FASTA File")
        #Elemen upload file
        uploaded_file = st.file_uploader("Drag FASTA file here!", type=["fasta", "fa"])
        
        if st.button("Run"):
            genome_data_with_header = uploaded_file.read().decode()

            #Deskripsi untuk line pertama  
            st.subheader("Description of the Genom")
            lines = genome_data_with_header.split("\n")
            if lines:
                desc = lines[0]
                st.write(desc)

            #Keseluruhan genom 
            st.subheader("Genome Data:")
            genome_data_no_header = "".join(genome_data_with_header.split("\n")[1:])  # Remove the first line (header)
            st.write(genome_data_no_header)


            st.subheader("Frequency each Nitrogen Base in The Sequence")
            Adenin = 0
            Guanin = 0
            Cytosin = 0
            Timin = 0 
            Error = 0
            for base in genome_data_no_header:
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
            
            st.write(f"Length of the data: {len(genome_data_no_header)}")
            data = {
                 "Parameter": ["Adenin", "Guanin", "Cytosin", "Timin", "Error rate"],
                 "Value (base)": [Adenin, Guanin, Cytosin, Timin, Error / len(genome_data_no_header)]
             }
            df = pd.DataFrame(data)
            st.table(df)
            
            st.subheader("Genome Metrics")

            # Create a pie chart for Adenin, Guanin, Timin, and Cytosin using Plotly
            fig = px.pie(
                names=['Adenin', 'Guanin', 'Timin', 'Cytosin'],
                values=[Adenin, Guanin, Timin, Cytosin],
                title='Genome Metrics: Adenin, Guanin, Timin, and Cytosin'
            )

            # Display the pie chart using Streamlit
            st.plotly_chart(fig)


if __name__ == "__main__":
    main()
