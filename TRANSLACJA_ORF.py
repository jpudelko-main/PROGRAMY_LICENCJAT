"""
Skrypt do translacji sekwencji nukleotydowych Salmonella i innych bakterii (katalogi A, B, C, D, R)
na sekwencje aminokwasowe, tworzący odpowiednią strukturę katalogów _aa z podkatalogami klasyfikacji.

Instrukcje:
1. Umieść ten skrypt w katalogu macierzystym, obok folderu 'SEKWENCJE_AA_KOMPLET' oraz pliku 'TABELA KODON NCBI.txt'.
2. Upewnij się, że masz prawa zapisu do tworzenia nowych katalogów i plików.
3. Uruchom skrypt: python3 translacja_spv.py
"""

import os
import sys

# --- KONFIGURACJA: ---
ROOT_DIR = "SEKWENCJE_AA_KOMPLET"           # Katalog źródłowy z podfolderami A, B, C, D, R
CODON_TABLE_FILE = "TABELA KODON NCBI.txt"   # Plik z tabelą kodonów NCBI

# Lista genów (nazwy katalogów w ROOT_DIR)
GENES = ["A", "B", "C", "D", "R"]

# Alternatywne kodony START (7)
START_CODONS = ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]
# Alternatywne kodony STOP (3)
STOP_CODONS = ["TAA", "TAG", "TGA"]

# Nazwy podkatalogów klasyfikacji w oryginale i odpowiadające im katalogi _aa
CLASS_ORIG = ["INNE_ORG", "SOLO_SEROWAR", "TOTAL"]
CLASS_AA   = ["INNE_ORG_aa", "SOLO_SEROWAR_aa", "TOTAL_aa"]

# -----------------------------------------------------------------------------
# Funkcja do wczytania i sparsowania tabeli kodonów NCBI:
def load_codon_table(path):
    """
    Wczytuje plik z tabelą kodonów NCBI w formacie:
      AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
      Starts = ---M------**--*----M------------MMMM---------------M------------
      Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
      Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
      Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

    Zwraca dwa słowniki:
      codon2aa: klucz = trzy nukleotydy (np. "ATG"), wartość = aminokwas (np. "M" lub "*")
      codon2start_flag: klucz = trójnukleotyd, wartość = True jeśli to kodon start (według linii Starts)
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]

    # Zakładamy, że plik ma co najmniej 5 niepuste linie
    # Pierwsza: "AAs    = ..."
    # Druga:   "Starts = ..."
    # Trzecia:  "Base1  = ..."
    # Czwarta: "Base2  = ..."
    # Piąta:   "Base3  = ..."
    aas_line    = None
    starts_line = None
    b1_line     = None
    b2_line     = None
    b3_line     = None

    for line in lines:
        if line.startswith("AAs"):
            aas_line = line.split("=", 1)[1].strip()
        elif line.startswith("Starts"):
            starts_line = line.split("=", 1)[1].strip()
        elif line.startswith("Base1"):
            b1_line = line.split("=", 1)[1].strip()
        elif line.startswith("Base2"):
            b2_line = line.split("=", 1)[1].strip()
        elif line.startswith("Base3"):
            b3_line = line.split("=", 1)[1].strip()

    if not (aas_line and starts_line and b1_line and b2_line and b3_line):
        print("Błąd: nie udało się wczytać wszystkich linii z pliku tabeli kodonów.", file=sys.stderr)
        sys.exit(1)

    if not (len(aas_line) == len(starts_line) == len(b1_line) == len(b2_line) == len(b3_line)):
        print("Błąd: długości wierszy tabeli kodonów nie są zgodne (musiałyby być 64).", file=sys.stderr)
        sys.exit(1)

    codon2aa = {}
    codon2start_flag = {}

    for i in range(len(aas_line)):
        codon = b1_line[i] + b2_line[i] + b3_line[i]
        aa    = aas_line[i]
        start_flag = (starts_line[i] == "M")  # True, jeśli w linii "Starts" jest "M"
        codon2aa[codon] = aa
        codon2start_flag[codon] = start_flag

    return codon2aa, codon2start_flag


# Funkcja czytająca sekwencję nukleotydową z pliku FASTA (przyjmuje ścieżkę do pliku .fasta)
def read_fasta_sequence(fasta_path):
    """
    Otwiera plik FASTA, ignoruje linię z nagłówkiem (zaczyna się od '>'),
    a następnie łączy wszystkie linie z sekwencją w jedną ciągłą sekwencję nukleotydową (uppercase).
    """
    seq_lines = []
    header = None
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line  # przechowujemy nagłówek (całą linię)
                # pomijamy ewentualne kolejne nagłówki, zakładamy jeden rekord na plik
                continue
            seq_lines.append(line.upper())
    sequence = "".join(seq_lines)
    return header, sequence


# Funkcja tłumacząca jedną sekwencję dla wybranej pary (start_codon, stop_codon),
# zwracająca sekwencję aminokwasową (bez znaku '*') lub pusty string, jeśli ORF niewystępuje
def translate_best_orf(nuc_seq, codon2aa, chosen_start, chosen_stop):
    """
    Dla danej sekwencji nukleotydowej:
    - Znajduje pierwsze wystąpienie chosen_start (TTG/CTG/...) w całym nuk_seq -> indeks i_start
    - Znajduje ostatnie wystąpienie chosen_stop (TAA/TAG/TGA) w nuk_seq -> indeks i_stop
      upewniając się, że i_stop > i_start. Jeśli takich kombinacji nie ma, zwraca pusty string.
    - ORF nukleotydowa = od i_start do i_stop (nie włącznie; stop nie jest tłumaczony)
    - Tłumaczy ORF od i_start, co 3 nukleotydy:
        * jeśli k == i_start: mapa AA = 'M'
        * w przeciwnym razie: mapa wg codon2aa[codon]
    - Ignoruje ewentualne niepełne końcowe nukleotydy, jeśli długość ORF nie jest wielokrotnością 3.
    - Zwraca string aminokwasowy (bez kodonu stop).
    """
    seq = nuc_seq
    i_start = seq.find(chosen_start)
    if i_start == -1:
        return ""

    i_stop = seq.rfind(chosen_stop)
    if i_stop == -1 or i_stop <= i_start:
        return ""

    # ORF nukleotydowy od i_start do i_stop (nie włącznie => nie tłumaczymy stopu)
    orf_nt = seq[i_start:i_stop]
    aa_seq = []

    # Długość ORF w nukleotydach
    L = len(orf_nt)
    # Tłumaczymy tylko całe triplety
    full_codons = (L // 3) * 3

    for pos in range(0, full_codons, 3):
        codon = orf_nt[pos:pos+3]
        if pos == 0:
            # Pierwszy kodon to zawsze chosen_start (sprawdziliśmy powyżej), wymuszamy 'M'
            aa_seq.append("M")
        else:
            aa = codon2aa.get(codon, "X")  # jeśli nieznany kodon, 'X'
            aa_seq.append(aa)

    return "".join(aa_seq)


# Główna funkcja przetwarzająca wszystkie pliki
def main():
    # 1. Wczytaj tabelę kodonów
    codon2aa, codon2start_flag = load_codon_table(CODON_TABLE_FILE)

    # 2. Sprawdź, czy istnieje katalog źródłowy
    if not os.path.isdir(ROOT_DIR):
        print(f"Błąd: nie znaleziono katalogu '{ROOT_DIR}'.", file=sys.stderr)
        sys.exit(1)

    # 3. Utwórz strukturę katalogów _aa: dla każdego genu A->A_aa, etc.
    for gene in GENES:
        src_gene_dir = os.path.join(ROOT_DIR, gene)
        if not os.path.isdir(src_gene_dir):
            print(f"Uwaga: nie znaleziono katalogu '{src_gene_dir}', pomijam ten gen.", file=sys.stderr)
            continue

        dst_gene_dir = os.path.join(ROOT_DIR, f"{gene}_aa")
        os.makedirs(dst_gene_dir, exist_ok=True)

        # Dla każdego podkatalogu klasyfikacji w źródle utwórz odpowiednik _aa
        for orig_cls, aa_cls in zip(CLASS_ORIG, CLASS_AA):
            src_cls_dir = os.path.join(src_gene_dir, orig_cls)
            dst_cls_dir = os.path.join(dst_gene_dir, aa_cls)
            if not os.path.isdir(src_cls_dir):
                print(f"Uwaga: nie znaleziono katalogu klasyfikacji '{src_cls_dir}', pomijam.", file=sys.stderr)
                continue
            os.makedirs(dst_cls_dir, exist_ok=True)

            # 4. Dla każdego pliku FASTA (.fasta) w src_cls_dir przetłumacz i zapisz w dst_cls_dir
            for fname in os.listdir(src_cls_dir):
                if not (fname.endswith(".fasta") or fname.endswith(".fa") or fname.endswith(".fna")):
                    continue

                src_fasta_path = os.path.join(src_cls_dir, fname)
                header, nuc_seq = read_fasta_sequence(src_fasta_path)
                if not nuc_seq:
                    print(f"Plik '{src_fasta_path}' jest pusty lub nieprawidłowy, pomijam.", file=sys.stderr)
                    continue

                best_aa = ""
                best_combo = None  # (start_codon, stop_codon)
                # 5. Próbuj wszystkich kombinacji start/stop
                for st in START_CODONS:
                    for sp in STOP_CODONS:
                        aa_seq = translate_best_orf(nuc_seq, codon2aa, st, sp)
                        if len(aa_seq) > len(best_aa):
                            best_aa = aa_seq
                            best_combo = (st, sp)

                # 6. Jeśli nie znaleziono ORF (best_aa == ""), to i tak zapisz pusty plik FASTA
                if best_combo is None:
                    # nie znaleziono żadnej kombinacji dającej ORF
                    best_combo = ("NA", "NA")  # opcjonalnie
                    best_aa = ""

                st_code, sp_code = best_combo
                # Usuń rozszerzenie z nazwą pliku
                name_root = os.path.splitext(fname)[0]
                # Nowa nazwa: oryginalna nazwa + start+stop + .fasta
                new_fname = f"{name_root}{st_code}{sp_code}.fasta"
                dst_fasta_path = os.path.join(dst_cls_dir, new_fname)

                # Zapisz FASTA
                with open(dst_fasta_path, "w", encoding="utf-8") as out_f:
                    if header:
                        out_f.write(header + "\n")
                    else:
                        out_f.write(f">{name_root}{st_code}{sp_code}\n")
                    # Zapis sekwencji: wiersze co 60 AA
                    for i in range(0, len(best_aa), 60):
                        out_f.write(best_aa[i:i+60] + "\n")

                print(f"Zapisano: {dst_fasta_path} (len AA = {len(best_aa)}, start={st_code}, stop={sp_code})")

    print("Przetwarzanie zakończone.")


if __name__ == "__main__":
    main()
