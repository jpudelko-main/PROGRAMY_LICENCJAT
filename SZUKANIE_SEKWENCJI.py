import os
import re
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez


Entrez.email = "j.pudelko.970@studms.ug.edu.pl"


# Definicja trzech faz i odpowiadających im zapytań Entrez
dphases = {
    "phase1": {"entrez_query": "Salmonella[Organism]"},
    "phase2": {"entrez_query": "Bacteria[Organism] NOT Salmonella[Organism]"},
    "phase3": {"entrez_query": "(Viruses[Organism] OR Fungi[Organism] OR Protista[Organism] OR Archaea[Organism] OR Plants[Organism])"}
}

# Pliki referencyjne: nukleotydowe i aminokwasowe spvr, spva, spvb,
# spvc, spvd - geny(lowercase) SpvR(p)... -> sekwencja aminokwasowa produktow

nt_reference_files = ["r.txt", "a.txt", "b.txt", "c.txt", "d.txt"]
aa_reference_files = ["Rp.txt", "Ap.txt", "Bp.txt", "Cp.txt", "Dp.txt"]

# Katalog bazowy na wyniki
base_output = "results"

# Minimalne pokrycie zapytania w procentach
def coverage_ok(hsp, query_len, threshold=80.0):
    return (hsp.align_length / query_len) * 100 >= threshold


def read_sequence(path):
    """
    Czyta plik .txt zawierający sekwencję referencyjną.
    Obsługuje zarówno FASTA (>...) jak i czysty tekst.
    """
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    if lines and lines[0].startswith('>'):
        return ''.join(lines[1:])
    return ''.join(lines)


def sanitize_name(name):
    """
    Czyści tytuł lub caption z niebezpiecznych znaków
    """
    cleaned = re.sub(r"[^\w\.-]", "_", name)
    return re.sub(r"_+", "_", cleaned).strip("_")


def process_results(handle, query_seq, phase, method, gene):
    seen = set()
    for rec in NCBIXML.parse(handle):
        for aln in rec.alignments:
            hsp = aln.hsps[0]
            if not coverage_ok(hsp, len(query_seq)):
                continue

            acc = aln.accession
            # Pobierz metadane przez Entrez
            with Entrez.esummary(db='nucleotide', id=acc) as sh:
                sm = Entrez.read(sh)[0]
            tid = sm.get('TaxId')
            if tid in seen:
                continue
            seen.add(tid)

            title = sm.get('Title') or sm.get('Caption') or acc
            name = sanitize_name(title)

            # Tworzenie katalogu wynikowego
            odir = os.path.join(base_output, phase, method, gene, f"{name}_{tid}")
            os.makedirs(odir, exist_ok=True)

            # Zapis pełnej sekwencji nukleotydowej
            with Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text') as fh:
                fasta = fh.read()
            with open(os.path.join(odir, 'full_sequence.fasta'), 'w') as f:
                f.write(fasta)
            so = ''.join(fasta.splitlines()[1:])
            with open(os.path.join(odir, 'full_sequence.txt'), 'w') as f:
                f.write(so)

            # Zapis fragmentu HSP
            sbj_start, sbj_end = min(hsp.sbjct_start, hsp.sbjct_end), max(hsp.sbjct_start, hsp.sbjct_end)
            frag = so[sbj_start-1:sbj_end]
            header = f">{acc}_{gene}_{method}_frag_{sbj_start}_{sbj_end}\n"
            with open(os.path.join(odir, 'fragment_sequence.fasta'), 'w') as f:
                f.write(header + frag)
            with open(os.path.join(odir, 'fragment_sequence.txt'), 'w') as f:
                f.write(frag)
            with open(os.path.join(odir, 'fragment_coords.txt'), 'w') as f:
                f.write(
                    f"Start: {sbj_start}\n"
                    f"End: {sbj_end}\n"
                    f"Coverage: {(hsp.align_length/len(query_seq))*100:.2f}%\n"
                )


def run_pipeline():
    for phase, info in dphases.items():
        query_filter = info['entrez_query']

        # Faza nukleotydowa: megablast i blastn
        for method in ['megablast', 'blastn']:
            mb_flag = (method == 'megablast')
            for ref in nt_reference_files:
                gene = os.path.splitext(ref)[0].upper()
                seq = read_sequence(ref)
                print(f"[{phase}] {method.upper()} dla genu {gene}")

                result_handle = NCBIWWW.qblast(
                    program='blastn',
                    database='nt',
                    sequence=seq,
                    entrez_query=query_filter,
                    hitlist_size=5000,
                    format_type='XML',
                    megablast=mb_flag
                )
                process_results(result_handle, seq, phase, method, gene)
                result_handle.close()

        # Faza aminokwasowa: tblastn
        method = 'tblastn'
        for ref in aa_reference_files:
            gene = os.path.splitext(ref)[0]
            seq = read_sequence(ref)
            print(f"[{phase}] {method.upper()} dla genu {gene}")

            result_handle = NCBIWWW.qblast(
                program='tblastn',
                database='nt',
                sequence=seq,
                entrez_query=query_filter,
                hitlist_size=5000,
                format_type='XML'
            )
            process_results(result_handle, seq, phase, method, gene)
            result_handle.close()


if __name__ == '__main__':
    run_pipeline()
