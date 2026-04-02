import logging
from collections import defaultdict

from .genome import GenomeProcessor
from .vcf import VCFProcessor
from .modifier import SequenceModifier 
from .protein import ProteinAnalyzer
from .utils import reverse_complement

class ResultsOutputter:
    VARIANT_SCHEMA = [
        "LOF_Info",
        "missense",
        "insertion",
        "deletion",
        "complex_indel",
        "exon_skip",
        "skipped_exons_detail",
        "intron_retention",
        "retained_introns_detail",
        "start_codon_loss",
        "stop_loss",
        "No_start_codon_for_reference",
        "same_as_other_transcript",
        "protein_loss",
        "ref_protein_empty",
        "alignment_failed"
    ]

    def __init__(self, genome_processor, vcf_processor, sequence_modifier, protein_analyzer):
        self.gp = genome_processor
        self.vp = vcf_processor
        self.sm = sequence_modifier
        self.pa = protein_analyzer
        self.header_written = False

    def _format_mutations(self, mutations):
        details = []
        for mut in mutations:
            labels_str = f"[{','.join(mut.labels)}]" if mut.labels else ""
            overlap_str = "(overlapping)" if mut.overlapping else ""
            details.append(f"{mut.id}{labels_str}{overlap_str}")
        return ";".join(details)

    def write_header(self, output_file, sample_file=None):
        if self.header_written:
            return

        header = [
           "Gene_ID", "Haplotype_ID", "mRNA", 
           "Haplotype_Count", "Frequency", "Variant_Type", 
           "Protein_Changes", "Haplotype_Mutations", 
           "Sample_Sources",  
           "Ref_Protein", "Alt_Protein", "Ref_CDS", "Alt_CDS",
           "Aligned_Ref", "Comparison_String", "Aligned_Alt" 
        ]
        output_file.write("\t".join(header) + "\n")
        
        if sample_file:
            sample_header = ["Gene_ID", "mRNA_ID", "Ref_ID", "Alt_IDs"] + self.vp.sample_names
            sample_file.write("\t".join(sample_header) + "\n")

        self.header_written = True

    def _build_fixed_variant_string(self, variant_types_set, structural_report):
        output_parts = []
        lof_str = ""
        same_as_str = ""
        for v in variant_types_set:
            if v.startswith("LOF(") or v.startswith("Non-LOF("):
                lof_str = v
            elif v.startswith("same_as_protein_of_"):
                same_as_str = v

        skipped_detail = ""
        if structural_report.get('skipped_exons'):
             items = [f"SkippedExon:[{x}]" for x in structural_report['skipped_exons']]
             skipped_detail = ";".join(items)
        retained_detail = ""
        if structural_report.get('retained_introns'):
             items = [f"RetainedIntron:[{x}]" for x in structural_report['retained_introns']]
             retained_detail = ";".join(items)

        for col_key in self.VARIANT_SCHEMA:
            val = ""            
            if col_key == "LOF_Info":
                val = lof_str
            elif col_key == "same_as_other_transcript":
                val = same_as_str
            elif col_key == "skipped_exons_detail":
                val = skipped_detail
            elif col_key == "retained_introns_detail":
                val = retained_detail
            else:
                if col_key in variant_types_set:
                    val = col_key
            output_parts.append(val)
        return "|".join(output_parts)

    def analyze_and_write_gene_results(self, gene_obj, gene_mrna_combinations, out_f, align_f, sample_f):
        if not gene_mrna_combinations:
            return

        gene_id = getattr(gene_obj, 'id', 'Unknown_Gene')
        total_haplotypes = len(self.vp.sample_names) * 2
        for mrna_id, combinations in gene_mrna_combinations.items():
            mrna_obj = self.gp.mRNA_data.get(mrna_id)
            if not mrna_obj: continue

            ref_cds = mrna_obj.cds_sequence
            if mrna_obj.strand == '-':
                ref_cds = reverse_complement(ref_cds)
            ref_protein, _ = self.pa.translate_cds(ref_cds, mrna_id)
            ref_hap_id = f"{mrna_id}:REF"
            unique_haplotypes_map = {} 
            unique_haplotypes_map[()] = {
                "id": ref_hap_id,
                "mutations": [],
                "index": 0,
                "count": 0,
                "samples": defaultdict(list)
            }
            sample_to_hap_indices = defaultdict(list) 
            local_hap_counter = 0

            for combo_id, combo_info in combinations.items():
                hap1_muts = combo_info['hap1_mutations']
                hap2_muts = combo_info['hap2_mutations']
                hap1_tuple = tuple(m.id for m in hap1_muts)
                hap2_tuple = tuple(m.id for m in hap2_muts)

                if hap1_tuple not in unique_haplotypes_map:
                    local_hap_counter += 1
                    unique_haplotypes_map[hap1_tuple] = {
                        "id": f"{mrna_id}:{local_hap_counter}",
                        "mutations": hap1_muts,
                        "index": local_hap_counter,
                        "count": 0,
                        "samples": defaultdict(list)
                    }
                
                if hap2_tuple not in unique_haplotypes_map:
                    local_hap_counter += 1
                    unique_haplotypes_map[hap2_tuple] = {
                        "id": f"{mrna_id}:{local_hap_counter}",
                        "mutations": hap2_muts,
                        "index": local_hap_counter,
                        "count": 0,
                        "samples": defaultdict(list)
                    }

                hap1_idx = unique_haplotypes_map[hap1_tuple]["index"]
                hap2_idx = unique_haplotypes_map[hap2_tuple]["index"]

                for sample_name, zygosity in combo_info['samples'].items():
                    unique_haplotypes_map[hap1_tuple]["count"] += 1
                    unique_haplotypes_map[hap2_tuple]["count"] += 1
                    
                    if zygosity == 'homo':
                        unique_haplotypes_map[hap1_tuple]["samples"][sample_name].append("Homo")
                        if sample_name not in sample_to_hap_indices:
                            sample_to_hap_indices[sample_name] = [hap1_idx, hap1_idx]
                    else:
                        unique_haplotypes_map[hap1_tuple]["samples"][sample_name].append("Hap1")
                        unique_haplotypes_map[hap2_tuple]["samples"][sample_name].append("Hap2")
                        if sample_name not in sample_to_hap_indices:
                            sample_to_hap_indices[sample_name] = [hap1_idx, hap2_idx]

            sorted_haps = sorted(unique_haplotypes_map.items(), key=lambda x: x[1]['index'])
            output_rows_buffer = []
            for hap_tuple, hap_info in sorted_haps:
                hap_muts_list = hap_info["mutations"]
                current_hap_id = hap_info["id"]
                count = hap_info["count"]                
                muts_str = self._format_mutations(hap_muts_list) if hap_muts_list else "."
                alt_cds, structural_report = self.sm.apply_variants_to_cds(mrna_obj, hap_muts_list)
                if mrna_obj.strand == '-':
                    alt_cds = reverse_complement(alt_cds)
                alt_protein, _ = self.pa.translate_cds(alt_cds, mrna_id)
                variant_types_set, changes, aligned_ref, aligned_alt = self.pa.analyze_protein_changes(
                    ref_protein, alt_protein, mrna_obj, hap_muts_list, structural_report, gene_obj
                )

                comparison_string = ""
                if aligned_ref and aligned_alt:
                    comp_chars = []
                    for r, a in zip(aligned_ref, aligned_alt):
                        if r == a:
                            comp_chars.append("|")
                        elif r != '-' and a != '-':
                            comp_chars.append("*") # Mismatch
                        else:
                            comp_chars.append(" ") # Indel gap
                    comparison_string = "".join(comp_chars)
                else:
                    if not aligned_ref: aligned_ref = ref_protein
                    if not aligned_alt: aligned_alt = alt_protein
                    if not comparison_string and aligned_ref == aligned_alt:
                        comparison_string = "|" * len(aligned_ref)
                
                var_type_str = self._build_fixed_variant_string(variant_types_set, structural_report)
                if "start_codon_loss" in variant_types_set:
                    alt_protein = ""

                if hap_info['index'] > 0:
                    align_f.write(f"Haplotype_ID: {current_hap_id}\n")
                    align_f.write(f"Gene: {gene_id} | mRNA: {mrna_id}\n")
                    align_f.write(f"Haplotype_Mutations: {muts_str}\n")
                    align_f.write(f"Variant_Type: {var_type_str}\n")
                    align_f.write(f"Protein_Changes: {changes}\n")
                    align_f.write(f"Alignment:\n")
                    align_f.write(f"  Ref: {aligned_ref}\n")
                    align_f.write(f"       {comparison_string}\n")
                    align_f.write(f"  Alt: {aligned_alt}\n\n")

                is_ref = (hap_info['index'] == 0)
                
                if is_ref:
                    display_count = "."
                    display_freq = "."
                    sample_sources_str = "."
                    total_samples_str = "."
                else:

                    frequency = count / total_haplotypes if total_haplotypes > 0 else 0
                    display_count = count
                    display_freq = f"{frequency:.6f}"

                    sample_details_list = []
                    all_samples_in_hap = []
                    for s_name, tags in hap_info['samples'].items():
                        all_samples_in_hap.append(s_name)
                        for tag in tags:
                            sample_details_list.append(f"{s_name}({tag})")
                        
                    sample_details_list.sort()
                    all_samples_in_hap.sort()
                    sample_sources_str = ";".join(sample_details_list)
                    total_samples_str = ";".join(all_samples_in_hap)

                row_data = {
                    "Gene_ID": gene_id,
                    "Haplotype_ID": current_hap_id,
                    "mRNA": mrna_id,
                    "Haplotype_Count": display_count,
                    "Frequency": display_freq,
                    "Variant_Type": var_type_str if var_type_str else "Reference",
                    "Protein_Changes": changes if changes else ".",
                    "Haplotype_Mutations": muts_str,
                    "Sample_Sources": sample_sources_str,
                    "Ref_Protein": ref_protein,
                    "Alt_Protein": alt_protein,
                    "Ref_CDS": ref_cds,
                    "Alt_CDS": alt_cds,
                    "Aligned_Ref": aligned_ref,
                    "Comparison_String": comparison_string,
                    "Aligned_Alt": aligned_alt
                }
                output_rows_buffer.append(row_data)

            header = [
               "Gene_ID", "Haplotype_ID", "mRNA", 
               "Haplotype_Count", "Frequency", "Variant_Type", 
               "Protein_Changes", "Haplotype_Mutations", 
               "Sample_Sources", 
               "Ref_Protein", "Alt_Protein", "Ref_CDS", "Alt_CDS",
               "Aligned_Ref", "Comparison_String", "Aligned_Alt" 
            ]
            
            for row in output_rows_buffer:
                line = [str(row.get(col, '')) for col in header]
                out_f.write("\t".join(line) + "\n")

            alt_ids = [h[1]["id"] for h in sorted_haps[1:]]
            alt_ids_str = ",".join(alt_ids) if alt_ids else "."
            
            matrix_row = [gene_id, mrna_id, ref_hap_id, alt_ids_str]
            for s_name in self.vp.sample_names:
                indices = sample_to_hap_indices.get(s_name)
                if indices:
                    gt_str = f"{indices[0]}|{indices[1]}"
                else:
                    gt_str = "0|0" 
                matrix_row.append(gt_str)
            sample_f.write("\t".join(matrix_row) + "\n")

