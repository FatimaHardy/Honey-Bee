import subprocess

subprocess.check_call(["pip", "install", "biopython"])

from Bio import Entrez

# Provide your email to NCBI (required)
Entrez.email = "your_email@example.com"

# Define the search term
search_term = "Apis mellifera [ORGN] AND mitochondrial AND complete genome"

# Search GenBank
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
record = Entrez.read(handle)
handle.close()

# Get the list of GenBank IDs
genbank_ids = record["IdList"]
print(f"Found {len(genbank_ids)} records:")
print(genbank_ids)

from Bio import SeqIO

# Fetch sequence data for the first GenBank ID
genbank_id = genbank_ids[0]  # You can loop over all IDs
handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
sequence_record = SeqIO.read(handle, "genbank")
handle.close()

# Print some details about the sequence
print(f"Sequence ID: {sequence_record.id}")
print(f"Sequence Description: {sequence_record.description}")
print(f"Sequence Length: {len(sequence_record.seq)}")

search_term = "Apis mellifera [ORGN] AND Europe [Location]"

SeqIO.write(sequence_record, "honeybee_sequence.fasta", "fasta")

from Bio import SeqIO

# Assume you have already retrieved a sequence record
# For example:
# handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
# sequence_record = SeqIO.read(handle, "genbank")

# Extract specific gene sequences
gene_sequences = []
for feature in sequence_record.features:
    if feature.type == "gene":
        gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
        gene_seq = feature.extract(sequence_record.seq)
        gene_sequences.append((gene_name, gene_seq))

# Print out the genes and their sequences
for gene_name, gene_seq in gene_sequences:
    print(f"Gene: {gene_name}\nSequence: {gene_seq[:50]}...")  # Print first 50 bases
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

# Write the gene sequences to a FASTA file for alignment
with open("genes.fasta", "w") as output_handle:
    for gene_name, gene_seq in gene_sequences:
        output_handle.write(f">{gene_name}\n{gene_seq}\n")

import subprocess
from Bio import SeqIO
from Bio import AlignIO

# Save your sequences in a FASTA file
with open("genes.fasta", "w") as output_handle:
    for gene_name, gene_seq in gene_sequences:
        output_handle.write(f">{gene_name}\n{gene_seq}\n")

# Call MAFFT from the command line using subprocess
subprocess.run(["mafft", "--auto", "genes.fasta"], stdout=open("aligned_genes.fasta", "w"))

# Read and print the alignment using Biopython's AlignIO
alignment = AlignIO.read("aligned_genes.fasta", "fasta")
print(alignment)

import allel
import numpy as np

# Load the VCF file
vcf_path = "example_data.vcf"
callset = allel.read_vcf(vcf_path)

# Extract genotype data (assuming diploid organisms)
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# Calculate allele counts
allele_counts = genotypes.count_alleles()

# Calculate allele frequencies
allele_freqs = allele_counts.to_frequencies()

print(allele_freqs)

import matplotlib.pyplot as plt

# Example: Plot allele frequencies for the first SNP
plt.hist(allele_freqs[:, 0], bins=np.arange(0, 1.1, 0.1))
plt.title("Allele Frequency Distribution for SNP 1")
plt.xlabel("Allele Frequency")
plt.ylabel("Count")
plt.show()

# Perform PCA
coords, model = allel.pca(genotypes.to_n_alt(), n_components=2)

# Plot the first two principal components
plt.scatter(coords[:, 0], coords[:, 1])
plt.title("PCA of Genetic Data")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()

