import logging
import re
import pysam
from collections import defaultdict, namedtuple
from Bio.Seq import Seq

from .utils import reverse_complement

Gene = namedtuple('Gene', ['id', 'chrom', 'start', 'end', 'strand', 'mRNAs'])
mRNA = namedtuple('mRNA', ['id', 'chrom', 'start', 'end', 'strand', 'regions', 'cds_sequence', 'splice_junctions'])
Region = namedtuple('Region', ['type', 'chrom', 'start', 'end', 'strand', 'sequence'])
Mutation = namedtuple('Mutation', ['id', 'chrom', 'pos', 'ref', 'alt', 'type', 'overlapping', 'labels'])

class GenomeProcessor:
    def __init__(self, gff_file, fasta_file):
        self.gff_file = gff_file
        self.fasta_file = fasta_file
        self.gene_data = {}
        self.mRNA_data = {}
        self.reference_sequences = {}
        self.chrom_lengths = {}

    def parse_gff(self):
        logging.info(f"Parsing GFF file: {self.gff_file}")
        gene_mRNAs = defaultdict(list)
        mrna_exons = defaultdict(list)
        mrna_cds = defaultdict(list)
        mrna_five_prime_utrs = defaultdict(list)
        mrna_three_prime_utrs = defaultdict(list)
        file_format = None
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                chrom = parts[0]
                feature_type = parts[2]
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                if '"' in parts[8] and ' ' in parts[8]:
                    if file_format is None:
                        logging.info("Detected GTF format.")
                        file_format = 'gtf'
                    attributes = dict(re.findall(r'(\w+)\s+"([^"]+)"', parts[8]))
                elif '=' in parts[8]:
                    if file_format is None:
                        logging.info("Detected GFF3 format.")
                        file_format = 'gff'
                    attributes = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
                else:
                    if file_format is None:
                        logging.warning("Could not definitively determine format. Assuming GFF3.")
                        file_format = 'gff'
                    attributes = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
                if file_format == 'gtf':
                    if feature_type == 'gene':
                        gene_id = attributes.get('gene_id')
                        if gene_id: self.gene_data[gene_id] = Gene(gene_id, chrom, start, end, strand, [])
                    elif feature_type == 'transcript' or feature_type == 'mRNA':
                        mrna_id, parent_gene_id = attributes.get('transcript_id'), attributes.get('gene_id')
                        if mrna_id and parent_gene_id:
                            self.mRNA_data[mrna_id] = mRNA(mrna_id, chrom, start, end, strand, defaultdict(list), '', [])
                            gene_mRNAs[parent_gene_id].append(mrna_id)
                    elif feature_type == 'exon':
                        parent_mrna_id = attributes.get('transcript_id')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_exons[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'CDS':
                        parent_mrna_id = attributes.get('transcript_id')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_cds[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'five_prime_UTR':
                        parent_mrna_id = attributes.get('transcript_id')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_five_prime_utrs[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'three_prime_UTR':
                        parent_mrna_id = attributes.get('transcript_id')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_three_prime_utrs[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                elif file_format == 'gff':
                    if feature_type == 'gene':
                        gene_id = attributes.get('ID')
                        if gene_id: self.gene_data[gene_id] = Gene(gene_id, chrom, start, end, strand, [])
                    elif feature_type == 'transcript' or feature_type == 'mRNA':
                        mrna_id, parent_gene_id = attributes.get('ID'), attributes.get('Parent')
                        if mrna_id and parent_gene_id:
                            self.mRNA_data[mrna_id] = mRNA(mrna_id, chrom, start, end, strand, defaultdict(list), '', [])
                            gene_mRNAs[parent_gene_id].append(mrna_id)
                    elif feature_type == 'exon':
                        parent_mrna_id = attributes.get('Parent')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_exons[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'CDS':
                        parent_mrna_id = attributes.get('Parent')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_cds[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'five_prime_UTR':
                        parent_mrna_id = attributes.get('Parent')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_five_prime_utrs[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
                    elif feature_type == 'three_prime_UTR':
                        parent_mrna_id = attributes.get('Parent')
                        if parent_mrna_id in self.mRNA_data:
                            mrna_three_prime_utrs[parent_mrna_id].append({'chrom': chrom, 'start': start, 'end': end})
        for mrna_id in mrna_exons: mrna_exons[mrna_id].sort(key=lambda x: x['start'])
        for mrna_id in mrna_cds: mrna_cds[mrna_id].sort(key=lambda x: x['start'])
        for mrna_id, mrna_obj in self.mRNA_data.items():
            for exon in mrna_exons.get(mrna_id, []):
                mrna_obj.regions['EXON'].append(Region('EXON', mrna_obj.chrom, exon['start'], exon['end'], mrna_obj.strand, ''))
            for cds in mrna_cds.get(mrna_id, []):
                 mrna_obj.regions['CDS'].append(Region('CDS', mrna_obj.chrom, cds['start'], cds['end'], mrna_obj.strand, ''))
            for utr in mrna_five_prime_utrs.get(mrna_id, []):
                mrna_obj.regions['five_prime_UTR'].append(Region('five_prime_UTR', chrom, utr['start'], utr['end'], strand, ''))
            for utr in mrna_three_prime_utrs.get(mrna_id, []):
                mrna_obj.regions['three_prime_UTR'].append(Region('three_prime_UTR', chrom, utr['start'], utr['end'], strand, ''))

            sorted_exons = sorted(mrna_obj.regions['EXON'], key=lambda r: r.start)
            if len(sorted_exons) > 1:
                for i in range(len(sorted_exons) - 1):
                    intron_start = sorted_exons[i].end + 1
                    intron_end = sorted_exons[i+1].start - 1
                    if intron_start <= intron_end:
                        mrna_obj.regions['INTRON'].append(Region('INTRON', mrna_obj.chrom, intron_start, intron_end, mrna_obj.strand, ''))

            sorted_introns = sorted(mrna_obj.regions['INTRON'], key=lambda r: r.start)
            filtered_introns = [i for i in sorted_introns if (i.end - i.start) >= 4]
            for intron in filtered_introns:
                chrom_len = self.chrom_lengths.get(mrna_obj.chrom, 0)
                if mrna_obj.strand == '+': # define donor and acceptor sites and their windows
                    donor_start, donor_end = intron.start, intron.start + 1
                    acceptor_start, acceptor_end = intron.end - 1, intron.end
                    donor_window_start = max(1, intron.start - 3)
                    donor_window_end = min(chrom_len, intron.start + 8)
                    acceptor_window_start = max(1, intron.end - 8)
                    acceptor_window_end = min(chrom_len, intron.end + 3)
                else:
                    donor_start, donor_end = intron.end - 1, intron.end
                    acceptor_start, acceptor_end = intron.start, intron.start + 1
                    donor_window_start = max(1, intron.end - 8)
                    donor_window_end = min(chrom_len, intron.end + 3)
                    acceptor_window_start = max(1, intron.start - 3)
                    acceptor_window_end = min(chrom_len, intron.start + 8)
                if donor_start <= donor_end:
                    donor_junction = {'type': 'donor', 'site': Region('splice_donor_site', mrna_obj.chrom, donor_start, donor_end, mrna_obj.strand, ''), 'window': Region('splice_donor_window', mrna_obj.chrom, donor_window_start, donor_window_end, mrna_obj.strand, '')}
                    mrna_obj.splice_junctions.append(donor_junction)
                else:
                    logging.warning(f"Invalid donor site coordinates for mRNA {mrna_obj.id}: start ({donor_start}) > end ({donor_end}). Skipping.")
                if acceptor_start <= acceptor_end:
                    acceptor_junction = {'type': 'acceptor', 'site': Region('splice_acceptor_site', mrna_obj.chrom, acceptor_start, acceptor_end, mrna_obj.strand, ''), 'window': Region('splice_acceptor_window', mrna_obj.chrom, acceptor_window_start, acceptor_window_end, mrna_obj.strand, '')}
                    mrna_obj.splice_junctions.append(acceptor_junction)
                else:
                    logging.warning(f"Invalid acceptor site coordinates for mRNA {mrna_obj.id}: start ({acceptor_start}) > end ({acceptor_end}). Skipping.")
        for gene_id, mrna_ids in gene_mRNAs.items():
            if gene_id in self.gene_data:
                self.gene_data[gene_id] = self.gene_data[gene_id]._replace(mRNAs=[self.mRNA_data[mid] for mid in mrna_ids if mid in self.mRNA_data])
        logging.info("GFF parsing complete.")

    def load_reference_sequences(self):
        logging.info(f"Loading reference sequences from FASTA file: {self.fasta_file}")
        try:
            with pysam.FastaFile(self.fasta_file) as fasta:
                for chrom_name in fasta.references:
                   self.reference_sequences[chrom_name] = fasta.fetch(chrom_name)
                   self.chrom_lengths[chrom_name] = fasta.get_reference_length(chrom_name)
            logging.info("Reference sequences loaded.")
        except Exception as e:
            logging.error(f"Error loading FASTA file: {e}")
            raise

    def extract_region_sequences(self):
        logging.info("Extracting sequences for GFF regions.")
        for mrna_obj in self.mRNA_data.values():
            for regions in mrna_obj.regions.values():
                for i, region in enumerate(regions):
                    if region.chrom not in self.reference_sequences: continue
                    seq = self.reference_sequences[region.chrom][region.start - 1:region.end]
                    regions[i] = region._replace(sequence=seq)
            if hasattr(mrna_obj, 'splice_junctions'):
                for i, junction in enumerate(mrna_obj.splice_junctions):
                    site_region = junction['site']
                    window_region = junction['window']
                    if site_region.chrom in self.reference_sequences:
                        site_seq = self.reference_sequences[site_region.chrom][site_region.start - 1:site_region.end]
                        mrna_obj.splice_junctions[i]['site'] = site_region._replace(sequence=site_seq)
                    if window_region.chrom in self.reference_sequences:
                        window_seq = self.reference_sequences[window_region.chrom][window_region.start - 1:window_region.end]
                        mrna_obj.splice_junctions[i]['window'] = window_region._replace(sequence=window_seq)
        logging.info("Region sequence extraction complete.")

    def assemble_cds_sequences(self):
        logging.info("Assembling CDS sequences for mRNAs.")
        for mrna_id, mrna_obj in self.mRNA_data.items():
            cds_regions = mrna_obj.regions['CDS']
            sorted_cds = sorted(cds_regions, key=lambda r: r.start)
            full_cds_seq = "".join(cds.sequence for cds in sorted_cds)
            self.mRNA_data[mrna_id] = mrna_obj._replace(cds_sequence=full_cds_seq)
        for gene_id, gene_obj in self.gene_data.items():
            updated_mRNAs = []
            for old_mrna in gene_obj.mRNAs:
                if old_mrna.id in self.mRNA_data:
                    updated_mRNAs.append(self.mRNA_data[old_mrna.id])
                else:
                    updated_mRNAs.append(old_mrna)
            new_gene_obj = gene_obj._replace(mRNAs=updated_mRNAs)
            self.gene_data[gene_id] = new_gene_obj
        logging.info("CDS sequence assembly complete.")

