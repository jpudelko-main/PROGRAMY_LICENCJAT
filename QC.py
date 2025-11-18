import os
import glob
import pandas as pd


def read_fasta_sequence(fasta_path):
    """
    Reads a FASTA file and returns the concatenated sequence (lines after the header).
    """
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    return ''.join(line.strip() for line in lines[1:])


def find_reference_length(ref_folder, gene):
    """
    Finds the reference FASTA file for a given gene in ref_folder and returns its length.
    Matches files containing the gene name.
    """
    pattern = os.path.join(ref_folder, f"*{gene}*.fasta")
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No reference FASTA found for gene '{gene}' in {ref_folder}")
    ref_path = files[0]
    print(f"Using reference file for gene '{gene}': {os.path.basename(ref_path)}")
    return len(read_fasta_sequence(ref_path))


def compare_gene_lengths(script_dir, output_file="gene_length_comparison.xlsx"):
    # Define folders
    ref_folder = os.path.join(script_dir, 'reff')
    cmp_folder = script_dir

    # Extract genes from reference folder filenames
    ref_files = glob.glob(os.path.join(ref_folder, '*.fasta'))
    genes = []
    for path in ref_files:
        name = os.path.basename(path)
        gene = name.split()[0]  # take first token as gene name
        genes.append(gene)

    results = {}
    for gene in genes:
        length_ref = find_reference_length(ref_folder, gene)
        results[gene] = []
        gene_folder = os.path.join(cmp_folder, gene)
        if not os.path.isdir(gene_folder):
            print(f"Warning: folder for gene '{gene}' not found: {gene_folder}")
            continue
        for fasta_file in os.listdir(gene_folder):
            if not fasta_file.lower().endswith('.fasta'):
                continue
            cmp_path = os.path.join(gene_folder, fasta_file)
            seq_cmp = read_fasta_sequence(cmp_path)
            percent = (len(seq_cmp) / length_ref) * 100
            organism = os.path.splitext(fasta_file)[0]
            results[gene].append((organism, percent))

    # Write results to Excel
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for gene, data in results.items():
            df = pd.DataFrame(data, columns=["Organism", f"{gene}_Length_Percent"])
            df.to_excel(writer, sheet_name=gene, index=False)

    print(f"Comparison complete. Results saved to '{output_file}'")


if __name__ == "__main__":
    # Determine script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    compare_gene_lengths(script_dir)
