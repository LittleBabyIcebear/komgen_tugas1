import streamlit as st
import pandas as pd
import random
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

    selected_option = st.sidebar.radio("Select an option", ["Show Amino Acid Table", "Analyze Sequence to Amino Acid"])

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

if __name__ == "__main__":
    main()
