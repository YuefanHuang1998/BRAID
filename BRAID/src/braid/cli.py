import sys
import argparse
import os
import pysam
import logging

from .utils import setup_logging, CheckpointManager
from .modifier import SequenceModifier 
from .protein import ProteinAnalyzer
from .genome import GenomeProcessor
from .vcf import VCFProcessor
from .output import ResultsOutputter
from .version import FancyVersionAction

try:
    from importlib import resources
except ImportError:
    import importlib_resources as resources

def get_test_data_path(filename):
    try:
        data_path = resources.files('braid').joinpath('data').joinpath(filename)
        return str(data_path)
    except AttributeError:
        with resources.path('braid', 'data') as data_dir:
            return str(data_dir / filename)

def run_test():
    vcf_path = get_test_data_path("test.vcf.gz")
    fasta_path = get_test_data_path("test.fasta")
    gff_path = get_test_data_path("test.gff3")
    test_args = [
        "-v", vcf_path,
        "-r", fasta_path,
        "-g", gff_path
    ]
    try:
        main(test_args)
    except Exception as e:
        print(f"\n Failed: {e}")
        import traceback
        traceback.print_exc()

def main(args_list=None):
    parser = argparse.ArgumentParser(
        description="Analyzes phased VCF data to predict variant effects on protein sequences for each haplotype.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-g','--gff', help='Input GFF3 file.')
    parser.add_argument('-r', '--reference', help='Reference genome FASTA file (indexed: .fai).')
    parser.add_argument('-v', '--vcf',  help='Phased VCF file (indexed: .tbi/.csi).')
    parser.add_argument('-o', '--output', default='variant_analysis_output.tsv', help='Output TSV file name.')
    parser.add_argument('--force-unphased', action='store_true', help='Skip the phasing check and force the script to run on a potentially unphased VCF file.')
    parser.add_argument('--ignore-intron', action='store_true', help='Ignore mutations marked only as intron.')
    parser.add_argument('-s', '--sample', help='Path to file with specific sample IDs to analyze (one per line, no header).')
    parser.add_argument('--gene', help='Path to file with specific gene IDs to analyze (one per line, no header).') 
    parser.add_argument('--lof-threshold', type=float, default=0.30, help="LOF classification threshold.")
    parser.add_argument('--resume', action='store_true', help='Resume from the last successfully processed gene found in the log.')
    parser.add_argument('--version', nargs=0, action=FancyVersionAction, help='Show advanced version and environment information.')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    parser_test = subparsers.add_parser("test", help="Run tests using built-in data.")

    args = parser.parse_args(args_list)

    if args.command == "test":
        run_test()
        return
    if not args.gff or not args.reference or not args.vcf:
        parser.print_help()
        logging.error("Error: Arguments -g, -r, and -v are required for analysis.")
        sys.exit(1)

    output_basename = os.path.splitext(args.output)[0]
    args.log = f"{output_basename}.log"
    setup_logging(args.log)
    
    args.align = f"{output_basename}.alignment.txt"
    args.sample_matrix = f"{output_basename}.sample.txt"

    ckpt_manager = CheckpointManager(args.log)
    open_mode = 'w'

    if args.resume:
        logging.info("Received (Resume Mode)...")
        ckpt_manager.load_completed_genes()
        ckpt_manager.truncate_if_incomplete(args.output, file_type='tsv')
        ckpt_manager.truncate_if_incomplete(args.sample_matrix, file_type='tsv')
        ckpt_manager.truncate_if_incomplete(args.align, file_type='align')
        
        open_mode = 'a'
    else:
        logging.info("Start new task (Start from scratch)...")

    try:
        genome = GenomeProcessor(args.gff, args.reference)
        genome.load_reference_sequences()
        genome.parse_gff()
        genome.extract_region_sequences()
        genome.assemble_cds_sequences()

        all_genes = list(genome.gene_data.items()) 
        genes_to_process = all_genes
        if args.gene:
            if not os.path.exists(args.gene):
                logging.error(f"Gene list file not found: {args.gene}")
                sys.exit(1)
            logging.info(f"Loading target genes from {args.gene}...")
            with open(args.gene, 'r') as f:
                target_gene_ids = set(line.strip() for line in f if line.strip())            
            genes_to_process = [(gid, obj) for gid, obj in all_genes if gid in target_gene_ids]            
            logging.info(f"Target genes loaded: {len(target_gene_ids)}. Found in GFF: {len(genes_to_process)}.")
            if len(genes_to_process) == 0:
                logging.warning("No target genes found in the GFF file! Please check IDs.")
                sys.exit(0)

        vcf_processor = VCFProcessor(genome, args.vcf, args.force_unphased, args.sample, args.ignore_intron)
        sequence_modifier = SequenceModifier(genome)
        protein_analyzer = ProteinAnalyzer(genome, args.lof_threshold)
        
        results_outputter = ResultsOutputter(genome, vcf_processor, sequence_modifier, protein_analyzer)

        if args.resume and os.path.exists(args.output) and os.path.getsize(args.output) > 0:
            results_outputter.header_written = True

        with pysam.VariantFile(args.vcf) as vcf_in, \
             open(args.output, open_mode, buffering=1) as out_f, \
             open(args.align, open_mode, buffering=1) as align_f, \
             open(args.sample_matrix, open_mode, buffering=1) as sample_f:

            vcf_processor.load_and_validate_samples(vcf_in)
            results_outputter.write_header(out_f, sample_f)
            logging.info(f"Starting pipeline. Total genes: {len(genes_to_process)}")

            for i, (gene_id, gene_obj) in enumerate(genes_to_process, 1):
                if gene_id in ckpt_manager.completed_genes:
                    continue
                logging.info(f"Processing Gene {i}/{len(genes_to_process)}: {gene_id}...")

                gene_variant_combinations = vcf_processor.process_variants_for_gene(gene_obj, vcf_in)
                if gene_variant_combinations:
                    results_outputter.analyze_and_write_gene_results(gene_obj, gene_variant_combinations, out_f, align_f, sample_f)

                    out_f.flush()
                    align_f.flush()
                    sample_f.flush()

                ckpt_manager.log_completion(gene_id)

        logging.info("Pipeline finished successfully.")

    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()

