import logging
import re
import os
import sys
import pysam
from collections import defaultdict, namedtuple

from .genome import mRNA, Mutation
from .utils import reverse_complement

class VCFProcessor:
    def __init__(self, genome_processor, vcf_file, force_unphased=False, sample_file=None, ignore_intron=False):
        self.genome_processor = genome_processor
        self.vcf_file = vcf_file
        self.force_unphased = force_unphased
        self.sample_names = []
        self.sample_file = sample_file
        self.logged_homo_ref_in_hap = set()
        self._vcf_samples_loaded = False
        self.ignore_intron = ignore_intron 

    def load_and_validate_samples(self, vcf_in):
        if self._vcf_samples_loaded:
            return
        logging.info("Loading and validating VCF samples...")

        vcf_path = vcf_in.filename
        if isinstance(vcf_path, bytes):
            vcf_path = vcf_path.decode('utf-8')
        if not os.path.exists(vcf_path + ".tbi") and not os.path.exists(vcf_path + ".csi"):
            error_msg = f"Index tbi/csi is missing, run 'bcftools index {vcf_path}' to produce index"
            print(error_msg)
            logging.error(error_msg)
            sys.exit(1)

        all_vcf_samples = list(vcf_in.header.samples)
        if self.sample_file:
            logging.info(f"Reading sample IDs from file: {self.sample_file}")
            with open(self.sample_file, 'r') as f:
                requested_samples = [line.strip() for line in f if line.strip()]
            vcf_samples_set = set(all_vcf_samples)
            missing_samples = [s for s in requested_samples if s not in vcf_samples_set]
            if missing_samples:
                raise ValueError(f"The following sample IDs from your file were not found in the VCF header: {', '.join(missing_samples)}")
            self.sample_names = requested_samples
            logging.info(f"Successfully loaded {len(self.sample_names)} samples for analysis.")
        else:
            self.sample_names = all_vcf_samples
            logging.info(f"No sample file provided. Analyzing all {len(self.sample_names)} samples found in VCF.")
        if not self.sample_names:
            raise ValueError("No samples selected for analysis. Please check your VCF file or sample list.")
        if not self.force_unphased:
            logging.info("Checking VCF for phased genotype information...")
            is_phased = True # Assume phased until proven otherwise, check first 1000 records
            records_to_check = 1000
            first_sample = self.sample_names[0]
            for i, record in enumerate(vcf_in):
                if i >= records_to_check: break
                if not record.samples[first_sample].phased:
                    is_phased = False
                    break
            vcf_in.reset()
            if not is_phased:
                error_msg = ("VCF file does not appear to be phased. This script requires phased genotypes (e.g., '0|1', not '0/1') to correctly construct haplotypes. Please provide a phased VCF.")
                logging.error(error_msg)
                raise ValueError(error_msg)
            logging.info("Phasing check passed. VCF appears to be phased.")
        self._vcf_samples_loaded = True

    def _label_variant(self, record, mrna_obj):
        labels = []
        if record.chrom != mrna_obj.chrom:
            logger.warning(f"Chromosome mismatch: VCF chrom {record.chrom} vs mRNA chrom {mrna_obj.chrom}")
            return []
        for region_type, regions in mrna_obj.regions.items():
            for region in regions:
                if record.pos <= region.end and record.stop >= region.start:
                    labels.append(region_type)
        if hasattr(mrna_obj, 'splice_junctions'):
            for junction in mrna_obj.splice_junctions:
                site_region = junction['site']
                if record.pos <= site_region.end and record.stop >= site_region.start:
                    labels.append(site_region.type)
        return sorted(list(set(labels)))

    def _create_mutation_from_record(self, record, alt_allele, labels):
        chrom = record.chrom
        ref = record.ref
        alt = '' if alt_allele is None or alt_allele == '*' else alt_allele
        mut_type = "SNP" if len(ref) == 1 and len(alt) == 1 else "INDEL"
        mutation_id = f"{chrom}:{record.pos}_{ref}>{alt}"
        return Mutation(mutation_id, chrom, record.pos, ref, alt, mut_type, False, labels)

    def _mark_overlaps(self, mutations, haplotype_id):
        if len(mutations) < 2: return
        max_end = mutations[0].pos + len(mutations[0].ref) - 1
        overlapping_indices = set()
        for i in range(1, len(mutations)):
            current_mut = mutations[i]
            prev_mut = mutations[i-1]
            if current_mut.pos <= max_end:
                overlapping_indices.add(i)
                overlapping_indices.add(i-1)
                short_hap_id = str(haplotype_id)[:50] + "..." if len(str(haplotype_id)) > 50 else str(haplotype_id) # Only log first 50 chars if too long
                logging.info(f"Overlapping mutations on same haplotype [[haplotype id: {short_hap_id}]]: {prev_mut.id} and {current_mut.id}")
            max_end = max(max_end, current_mut.pos + len(current_mut.ref) - 1)
        for i in overlapping_indices:
            if not mutations[i].overlapping:
                mutations[i] = mutations[i]._replace(overlapping=True)

    def process_variants_for_gene(self, gene_obj, vcf_in):
        gene_mrna_combinations = defaultdict(dict)
        sample_mrna_haplotypes = defaultdict(lambda: defaultdict(lambda: {'hap1': {}, 'hap2': {}}))

        try:
            fetched_records = list(vcf_in.fetch(gene_obj.chrom, gene_obj.start, gene_obj.end))
        except ValueError:
            logging.warning(f"Chromosome '{gene_obj.chrom}' for gene '{gene_obj.id}' not in VCF. Skipping.")
            return {}

        for record in fetched_records:
            pos = record.pos
            alleles = [record.ref] + list(record.alts)
            for mrna_obj in gene_obj.mRNAs:
                mrna_id = mrna_obj.id
                if not (record.pos <= mrna_obj.end and record.stop >= mrna_obj.start):
                    continue
                labels = self._label_variant(record, mrna_obj)
                if not labels: continue
                if self.ignore_intron and labels == ['INTRON']:
                    continue
                for sample_name in self.sample_names:
                    genotype = record.samples[sample_name].get('GT')
                    if genotype is None or len(genotype) < 2 or None in genotype: continue
                    hap_keys = ['hap1', 'hap2'] # Assuming diploid
                    for i in range(2):
                        allele_idx = genotype[i]
                        if allele_idx > 0 and alleles[allele_idx] is not None and alleles[allele_idx] != '*':
                            mut_obj = self._create_mutation_from_record(record, alleles[allele_idx], labels)
                            hap_key = hap_keys[i]
                            if pos not in sample_mrna_haplotypes[mrna_id][sample_name][hap_key]:
                                sample_mrna_haplotypes[mrna_id][sample_name][hap_key][pos] = mut_obj

        unique_haplotypes_for_mrna = defaultdict(dict)
        for mrna_id, samples_data in sample_mrna_haplotypes.items():
            for sample_name, haplotypes in samples_data.items():
                hap1_muts_list = sorted(haplotypes['hap1'].values(), key=lambda m: m.pos)
                hap2_muts_list = sorted(haplotypes['hap2'].values(), key=lambda m: m.pos)
                hap1_tuple = tuple(m.id for m in hap1_muts_list)
                hap2_tuple = tuple(m.id for m in hap2_muts_list)
                if hap1_tuple and hap1_tuple not in unique_haplotypes_for_mrna[mrna_id]:
                    unique_haplotypes_for_mrna[mrna_id][hap1_tuple] = hap1_muts_list
                if hap2_tuple and hap2_tuple not in unique_haplotypes_for_mrna[mrna_id]:
                    unique_haplotypes_for_mrna[mrna_id][hap2_tuple] = hap2_muts_list

                combo_representation = tuple((hap1_tuple, hap2_tuple))
                combo_id = str(combo_representation)
                if combo_id not in gene_mrna_combinations[mrna_id]:
                    gene_mrna_combinations[mrna_id][combo_id] = {
                        'hap1_tuple': hap1_tuple,
                        'hap2_tuple': hap2_tuple,
                        'samples': defaultdict(str)
                    }
                zygosity = 'homo' if hap1_tuple == hap2_tuple else 'hetero'
                gene_mrna_combinations[mrna_id][combo_id]['samples'][sample_name] = zygosity
        # Identify overlaps within haplotype
        for mrna_id, unique_haps in unique_haplotypes_for_mrna.items():
            for hap_tuple, muts_list in unique_haps.items():
                self._mark_overlaps(muts_list, hap_tuple)

        for mrna_id, combos in gene_mrna_combinations.items():
            for combo_id, combo_data in combos.items():
                hap1_tuple = combo_data.pop('hap1_tuple')
                hap2_tuple = combo_data.pop('hap2_tuple')

                combo_data['hap1_mutations'] = unique_haplotypes_for_mrna[mrna_id].get(hap1_tuple, [])
                combo_data['hap2_mutations'] = unique_haplotypes_for_mrna[mrna_id].get(hap2_tuple, [])
        
        return gene_mrna_combinations

