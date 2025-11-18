

import os
import sys
import time
import requests
import xml.etree.ElementTree as ET

# --- KONFIGURACJA ---
FOLDER   = "SOLO SEROWAR"
EMAIL    = "j.pudelko.970@studms.ug.edu.pl"
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"

# --- 1. Wczytywanie sekwencji FASTA ---
def load_fasta_sequences(folder):
    seqs = []
    for fname in os.listdir(folder):
        if fname.lower().endswith((".fasta", ".fa", ".fas")):
            path = os.path.join(folder, fname)
            with open(path, "r", encoding="utf-8") as f:
                seqs.append(f.read().strip())
    return "\n".join(seqs)

# --- 2. Wysyłanie zadania ---
def submit_job(fasta_data):
    r = requests.post(f"{BASE_URL}/run/", data={
        "email": EMAIL,
        "sequence": fasta_data
    })
    r.raise_for_status()
    return r.text.strip()

# --- 3. Polling statusu ---
def wait_for_completion(job_id, interval=5):
    while True:
        r = requests.get(f"{BASE_URL}/status/{job_id}")
        r.raise_for_status()
        status = r.text.strip()
        print(f"[{time.strftime('%H:%M:%S')}] Status: {status}")
        if status == "FINISHED":
            return
        if status not in ("RUNNING", "QUEUED"):
            raise RuntimeError(f"Zadanie zakończyło się błędem: {status}")
        time.sleep(interval)

# --- 4. Parsowanie dostępnych typów wyników ---
def get_result_types(job_id):
    r = requests.get(f"{BASE_URL}/resulttypes/{job_id}")
    r.raise_for_status()
    root = ET.fromstring(r.text)
    types = []
    for t in root.findall(".//type"):
        ident = t.findtext("identifier", default="").strip()
        desc  = t.findtext("description", default="").strip()
        if ident:
            types.append((ident, desc))
    return types

# --- 5. Pobranie i zapis wyniku ---
def fetch_result(job_id, result_type, out_path):
    r = requests.get(f"{BASE_URL}/result/{job_id}/{result_type}")
    if r.status_code == 400:
        raise ValueError(f"Niepoprawny format wyniku: {result_type}")
    r.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(r.content)
    print(f"✔ Wynik zapisano w: {out_path}")

# === GŁÓWNY BLOK ===
if __name__ == "__main__":
    fasta_data = load_fasta_sequences(FOLDER)
    if not fasta_data:
        print(f"Brak plików FASTA w folderze '{FOLDER}'.")
        sys.exit(1)

    print("1) Wysyłanie zadania…")
    job_id = submit_job(fasta_data)
    print("   Job ID:", job_id)

    print("2) Oczekiwanie na wynik…")
    wait_for_completion(job_id)

    print("3) Pobieram listę formatów wyników…")
    types = get_result_types(job_id)
    for ident, desc in types:
        print(f"  • {ident}: {desc}")

    # --- Automatyczny wybór formatu ZIP ---
    zip_id = None
    for ident, _ in types:
        if ident.lower() == "zip":
            zip_id = ident
            break

    if not zip_id:
        print("❌ Nie znaleziono formatu ZIP w dostępnych wynikach.")
        sys.exit(1)

    out_file = os.path.join(FOLDER, "msa_results.zip")
    print(f"4) Pobieram pakiet wyników ('{zip_id}')…")
    try:
        fetch_result(job_id, zip_id, out_file)
    except ValueError as e:
        print("❌", e)
        sys.exit(1)
