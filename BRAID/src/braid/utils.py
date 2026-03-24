import logging

def setup_logging(log_file="annotation.log"):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file), 
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging configured. Log file: {log_file}")

# Sequence reverse and genetic code
def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

class CheckpointManager:
    def __init__(self, log_file):
        self.log_file = log_file
        self.completed_genes = set()

    def load_completed_genes(self):
        if not os.path.exists(self.log_file):
            return
        with open(self.log_file, 'r') as f:
            for line in f:
                if "Finished processing gene:" in line:
                    parts = line.strip().split("Finished processing gene:")
                    if len(parts) > 1:
                        gene_id = parts[1].strip()
                        self.completed_genes.add(gene_id)
        logging.info(f"[CheckpointManager] already loaded {len(self.completed_genes)} completed genes from log.")

    def truncate_if_incomplete(self, file_path, file_type='tsv'):
        if not os.path.exists(file_path):
            return
        logging.info(f"[CheckpointManager] is checking file: {os.path.basename(file_path)}...")
        
        with open(file_path, 'r+') as f:
            header_pos = 0
            line = f.readline()
            if not line: return 
            header_pos = f.tell()

            while True:
                current_pos = f.tell() 
                line = f.readline()
                if not line: 
                    break
                
                current_gene_id = None
                try:
                    if file_type == 'tsv':
                        parts = line.split('\t')
                        if parts:
                            current_gene_id = parts[0]
                    elif file_type == 'align':
                        if line.startswith("Gene:"):
                            current_gene_id = line.split("|")[0].split(":")[1].strip()
                except Exception:
                    pass

                if current_gene_id and current_gene_id not in self.completed_genes:
                    logging.info(f"Find uncompleted data (Gene: {current_gene_id})，currently truncating file...")
                    f.seek(current_pos)
                    f.truncate()
                    break
        logging.info(f"Check completed: {os.path.basename(file_path)}")

    def log_completion(self, gene_id):
        logging.info(f"Finished processing gene: {gene_id}")

