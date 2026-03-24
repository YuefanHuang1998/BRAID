import logging
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from .genome import Gene, mRNA
from .utils import GENETIC_CODE

class ProteinAnalyzer:
    def __init__(self, genome_processor_instance, lof_threshold):
        self.genome_processor = genome_processor_instance
        self.aligner = PairwiseAligner()
        self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.mode = 'global'
        self.lof_threshold = lof_threshold
        logging.info("ProteinAnalyzer initialized with full, integrated functionality.")

    def translate_cds(self, cds_sequence, mrna_id):
        protein_sequence = []
        warning_messages = []
        
        if not cds_sequence:
            return "", ["Empty CDS sequence."]

        start_codon = cds_sequence[0:3].upper()
        if start_codon != "ATG":
            msg = f"Warning: Non-canonical start codon '{start_codon}' for mRNA {mrna_id}."
            warning_messages.append(msg)
            logging.warning(msg)
        
        for i in range(0, len(cds_sequence) - len(cds_sequence) % 3, 3):
            codon = cds_sequence[i:i+3].upper()
            amino_acid = GENETIC_CODE.get(codon, 'X') 
            protein_sequence.append(amino_acid)
            if amino_acid == '*': 
                break
        
        return "".join(protein_sequence), warning_messages

    def analyze_protein_changes(self, ref_protein, alt_protein, mrna_obj, mutations, structural_changes, gene_obj):
        variant_types = set()
        aligned_ref, aligned_alt = "", ""

        if structural_changes.get('skipped_exons'):
            variant_types.add("exon_skip")
        if structural_changes.get('retained_introns'):
            variant_types.add("intron_retention")

        if not ref_protein or not alt_protein:
            variant_types.add("ref_protein_empty" if not ref_protein else "protein_loss")
            return variant_types, "", aligned_ref, aligned_alt

        if ref_protein == alt_protein:
            is_splice_mut = any('splice' in label for mut in mutations for label in mut.labels)
            aligned_ref = ref_protein
            aligned_alt = alt_protein

        alignments = self.aligner.align(Seq(ref_protein), Seq(alt_protein))
        try:
            best_alignment = next(alignments)
            aligned_ref, aligned_alt = str(best_alignment[0]), str(best_alignment[1])
        except StopIteration:
            variant_types.add("alignment_failed")
            return variant_types, "", aligned_ref, aligned_alt

        changes = []
        ref_idx = 0
        i = 0
        identical_count = 0 

        while i < len(aligned_ref):
            ref_char, alt_char = aligned_ref[i], aligned_alt[i]

            if ref_char == alt_char and ref_char != '-':
                identical_count += 1

            if ref_char == '-' or alt_char == '-':
                indel_start_pos = ref_idx + 1
                ref_indel_seq, alt_indel_seq = "", ""
                
                while i < len(aligned_ref) and (aligned_ref[i] == '-' or aligned_alt[i] == '-'):
                    if aligned_ref[i] != '-': ref_indel_seq += aligned_ref[i]
                    if aligned_alt[i] != '-': alt_indel_seq += aligned_alt[i]
                    i += 1

                if ref_indel_seq and not alt_indel_seq: 
                    changes.append(f"Del({indel_start_pos}){ref_indel_seq}")
                    variant_types.add("deletion")
                elif not ref_indel_seq and alt_indel_seq:
                    changes.append(f"Ins({indel_start_pos - 1}){alt_indel_seq}")
                    variant_types.add("insertion")
                else:
                    changes.append(f"Complex({indel_start_pos}){ref_indel_seq}>{alt_indel_seq}")
                    variant_types.add("complex_indel")
                ref_idx += len(ref_indel_seq)
                continue

            else: 
                if ref_idx == 0 and ref_char != 'M':
                   variant_types.add("No_start_codon_for_reference")
                if ref_char != alt_char:
                    changes.append(f"Sub({ref_idx + 1}){ref_char}>{alt_char}")
                    if ref_idx == 0 and ref_char == 'M':
                        variant_types.add("start_codon_loss")
                    elif ref_char == '*': variant_types.add("stop_loss")
                    else: variant_types.add("missense")
                ref_idx += 1
                i += 1

        if ref_protein:
            ref_protein_len = len(ref_protein)
            non_identical_AAs = ref_protein_len - identical_count

            if "start_codon_loss" in variant_types: 
                non_identical_AAs = ref_protein_len
    
            non_identity_percentage = (non_identical_AAs / ref_protein_len) if ref_protein_len > 0 else 0
        
            if non_identity_percentage >= self.lof_threshold:
                lof_info = (f"LOF(non_identity_rate:{non_identity_percentage:.2%},"
                            f"non_identical_AAs:{non_identical_AAs},"
                            f"total_ref_AAs:{ref_protein_len})")
            else:
                lof_info = (f"Non-LOF(non_identity_rate:{non_identity_percentage:.2%},"
                            f"non_identical_AAs:{non_identical_AAs},"
                            f"total_ref_AAs:{ref_protein_len})")
            variant_types.add(lof_info)

        if gene_obj:
            for other_mrna_obj in gene_obj.mRNAs:
                if other_mrna_obj.id == mrna_obj.id:
                    continue

                other_ref_protein, _ = self.translate_cds(other_mrna_obj.cds_sequence, other_mrna_obj.id)
                if alt_protein and alt_protein == other_ref_protein:
                    variant_types.add(f"same_as_protein_of_{other_mrna_obj.id}")
                    logging.info(f"Protein from {mrna_obj.id} with variants is identical to reference protein of {other_mrna_obj.id}.")
        
        return variant_types, ";".join(changes), aligned_ref, aligned_alt

