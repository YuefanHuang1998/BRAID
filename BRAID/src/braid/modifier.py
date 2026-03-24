import logging
from collections import defaultdict
from Bio.Align import PairwiseAligner

from .genome import Gene, mRNA, Region
from .utils import reverse_complement

class SequenceModifier:
    def __init__(self, genome_processor):
        self.genome_processor = genome_processor
        self.logged_splice_events = set()
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5 

        self.aligner.target_left_open_gap_score = -10
        self.aligner.target_left_extend_gap_score = -10
        self.aligner.query_left_open_gap_score = -10
        self.aligner.query_left_extend_gap_score = -10

        self.aligner.target_right_open_gap_score = 0
        self.aligner.target_right_extend_gap_score = 0
        self.aligner.query_right_open_gap_score = 0
        self.aligner.query_right_extend_gap_score = 0

    def _get_alt_slice_for_region_with_alignment(self, ref, alt, var_start, region_start, region_end):
        var_end = var_start + len(ref)
        if var_end <= region_start or var_start >= region_end:
            return ""
        if var_start >= region_start and var_end <= region_end:
            return alt

        if not ref or not alt:
            logging.error("Reference or alternate sequence is empty. Skipping alignment.")
            return ""

        alignments = self.aligner.align(ref, alt)
        if not alignments:
            logging.error(f"Can't compare REF '{ref}' and ALT '{alt}' and position {var_start}")
            return ""

        best_alignment = alignments[0]
        aligned_ref, aligned_alt = str(best_alignment[0]), str(best_alignment[1])

        ref_bases_to_walk_before_region = max(0, region_start - var_start)
        ref_bases_to_walk_in_region = max(0, min(var_end, region_end) - max(var_start, region_start) + 1)

        ref_walk_count = 0
        alt_walk_count = 0
        alt_start_index = -1
        alt_end_index = -1

        for i in range(len(aligned_ref)):
            if ref_walk_count == ref_bases_to_walk_before_region and alt_start_index == -1:
                alt_start_index = alt_walk_count
            if aligned_ref[i] != '-':
                ref_walk_count += 1
            if aligned_alt[i] != '-':
                alt_walk_count += 1
            if ref_walk_count >= ref_bases_to_walk_before_region + ref_bases_to_walk_in_region:
                alt_end_index = alt_walk_count
                break

        if alt_end_index == -1:
            alt_end_index = len(alt)
        if alt_start_index == -1:
             alt_start_index = 0

        return alt[alt_start_index:alt_end_index]

    def analyze_splice_site_impact(self, mrna_obj, mutation, splice_window_region, original_site_region, mrna_strand):
        original_dinucleotide = original_site_region.sequence.upper() # Get the ACTUAL original dinucleotide from the 2bp site region
        site_len = len(original_dinucleotide)
        ref_window_seq = splice_window_region.sequence
        win_start = splice_window_region.start
        win_end = splice_window_region.start + len(ref_window_seq)

        site_start = original_site_region.start
        site_end = site_start + site_len

        mut_start = mutation.pos
        mut_ref = mutation.ref
        mut_alt = mutation.alt
        mut_end = mut_start + len(mut_ref)

        overlap_start_in_window = max(0, mut_start - win_start)
        overlap_end_in_window = min(len(ref_window_seq), mut_end - win_start)
        ref_seq_from_window = ref_window_seq[overlap_start_in_window:overlap_end_in_window]

        overlap_start_in_mut_ref = max(0, win_start - mut_start)
        overlap_end_in_mut_ref = min(len(mut_ref), win_end - mut_start)
        ref_seq_from_mutation = mut_ref[overlap_start_in_mut_ref:overlap_end_in_mut_ref]

        if ref_seq_from_window != ref_seq_from_mutation:
            logging.warning(
                f"Reference mismatch for mutation {mutation.id}. "
                f"Window slice '{ref_seq_from_window}' does not match "
                f"mutation.ref slice '{ref_seq_from_mutation}'. Skipping splice site analysis."
            )
            return 'REF_MISMATCH', None 

        prefix = ref_window_seq[:max(0, mut_start - win_start)]
        suffix = ref_window_seq[min(len(ref_window_seq), max(0, mut_end - win_start)):]
        mut_alt_for_window = self._get_alt_slice_for_region_with_alignment(
            ref=mut_ref, alt=mut_alt, var_start=mut_start,
            region_start=win_start, region_end=win_start + len(ref_window_seq)
        )
        alt_window_seq = prefix + mut_alt_for_window + suffix

        if not alt_window_seq:
            impact = 'SITE_DESTROYED'
            logging.info(
                f"[{mrna_obj.id}] (strand {mrna_obj.strand}) Splice window sequence became empty "
                f"due to mutation {mutation.id}. "
                f"Original site '{original_dinucleotide}' at {site_start} is DESTROYED."
            )
            return impact, site_start

        alt_dinucleotide = ""
        alignments = self.aligner.align(ref_window_seq, alt_window_seq)
        if not alignments:
            logging.error(f"Could not align window sequences for mutation {mutation.id}")
            return "ALIGNMENT_ERROR", None

        aligned_ref_win, aligned_alt_win = str(alignments[0][0]), str(alignments[0][1])
        local_site_start_in_ref = site_start - win_start
        
        align_end_idx = local_site_start_in_ref + site_len
        alt_dinucleotide = alt_window_seq[local_site_start_in_ref:align_end_idx]

        if alt_dinucleotide == original_dinucleotide:
            impact = 'SITE_PRESERVED'
            logging.info(f"[{mrna_obj.id}] (strand {mrna_obj.strand}) Splice site '{original_dinucleotide}' at {site_start} PRESERVED by mutation {mutation.id}.\n"
                         f"  - Original Window Seq  : {ref_window_seq}     - Mutated Window Seq  : {alt_window_seq}.")
            return impact, site_start

        new_found_pos = alt_window_seq.upper().find(original_dinucleotide)

        if new_found_pos != -1:
            impact = 'SITE_SHIFT'
            logging.warning(
                f"[{mrna_obj.id}] Splice site may SHIFT for mutation {mutation.id}.\n"
                f"  - Splice Site          : '{original_dinucleotide}' (Original genomic pos: {site_start})\n"
                f"  - Mutation             : {mutation.ref} -> {mutation.alt} (Genomic pos: {mutation.pos})\n"
                f"  - Window               : {splice_window_region.chrom}:{splice_window_region.start}-{splice_window_region.end}\n"
                f"  - Original Window Seq  : {ref_window_seq}\n"
                f"  - Mutated Window Seq   : {alt_window_seq}\n")
            return impact, site_start

        impact = 'SITE_DESTROYED'
        logging.info(
            f"[{mrna_obj.id}] (strand {mrna_obj.strand}) Splice site '{original_dinucleotide}' at {site_start} DESTROYED by mutation {mutation.id}. It was '{original_dinucleotide}', became '{alt_dinucleotide}'.\n"
            f"  - Original Window Seq  : {ref_window_seq}     - Mutated Window Seq  : {alt_window_seq}")
        return impact, site_start

    def get_alt_slice_for_cds(self, ref, alt, variant_start, region_start, region_end):
        if not ref or not alt:
            logging.error("Reference or alternate sequence is empty. Skipping alignment.")
            return ""

        variant_end = variant_start + len(ref)

        if variant_end <= region_start or variant_start > region_end:
            return ""  
        if variant_start >= region_start and variant_end <= (region_end + 1):
            return alt 

        alignments = self.aligner.align(ref, alt)
        if not alignments:
            logging.error(f"Could not align REF='{ref}' and ALT='{alt}' at pos {variant_start}")
            return "" 

        aligned_ref, aligned_alt = alignments[0]

        ref_bases_before_region = max(0, region_start - variant_start)
        ref_bases_in_region = max(0, min(variant_end, region_end) - max(variant_start, region_start) + 1)

        ref_walk_count = 0
        alt_walk_count = 0
        alt_start_index = -1
        alt_end_index = -1

        for i in range(len(aligned_ref)):
            if ref_walk_count >= ref_bases_before_region and alt_start_index == -1:
                alt_start_index = alt_walk_count
            if ref_walk_count >= ref_bases_before_region + ref_bases_in_region:
                alt_end_index = alt_walk_count
                break

            if aligned_ref[i] != '-':
                ref_walk_count += 1
            if aligned_alt[i] != '-':
                alt_walk_count += 1

        if alt_start_index == -1:
            alt_start_index = 0
        if alt_end_index == -1:
            alt_end_index = len(alt)

        logging.info(f"alt_start_index: {alt_start_index}, alt_end_index: {alt_end_index}")
        untruncated_alt_for_cds  = alt[alt_start_index:alt_end_index]

        return untruncated_alt_for_cds 

    def apply_variants_to_cds(self, mrna_obj, mutations):
        structural_report = {'skipped_exons': [], 'retained_introns': [], 'try_to_skip_exons': [], 'try_to_retain_introns': []}
        active_parts = sorted(mrna_obj.regions['CDS'], key=lambda r: r.start)        
        splice_mutations = [m for m in mutations if 'splice_donor_site' in m.labels or 'splice_acceptor_site' in m.labels]

        for mut in splice_mutations:

            splice_type = set()
            if 'splice_donor_site' in mut.labels:
                splice_type.add('donor')
            if 'splice_acceptor_site' in mut.labels: 
                splice_type.add('acceptor')

            impact_events = []
            for junction in mrna_obj.splice_junctions:
                mut_end = mut.pos + len(mut.ref) - 1
                window = junction['window']
                is_overlapping = not (mut_end < window.start or mut.pos > window.end)

                if junction['type'] in splice_type and is_overlapping:
                    target_window = junction['window']
                    target_site = junction['site']
                    impact, site_start = self.analyze_splice_site_impact(mrna_obj, mut, target_window, target_site, mrna_obj.strand)
                    if site_start is None:
                        continue

                    impact_events.append({
                        'impact': impact,
                        'site_start': site_start,
                        'junction_type': junction['type'],
                        'window': target_window,
                        'site': target_site,
                        'chrom': target_window.chrom
                        })
            
                    logging.info(f"Recorded impact '{impact}' for {junction['type']} at {site_start}")

            for event in impact_events:
                impact = event['impact']
                site_start = event['site_start']
                site_type = event['junction_type']                    
                if impact == 'SITE_PRESERVED' or site_start is None:
                    logging.info(f"Skipping structural change for mutation {mut.id} as its impact is '{impact}'.")
                    continue

                if site_type == 'acceptor':
                    logging.info(f"Try to simulating {mut}")
                    exon_to_skip = None
                    for exon in mrna_obj.regions['EXON']:
                        if mrna_obj.strand == '+' and exon.start - 2 <= site_start <= exon.start + 2:
                            exon_to_skip = exon
                            structural_report['try_to_skip_exons'].append(f"{exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end}|{impact}")
                            logging.info(f"add {exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end} to structural_report try_to_skip_exons")
                            break
                        elif mrna_obj.strand == '-' and exon.end - 2 <= site_start <= exon.end + 2:
                            exon_to_skip = exon
                            structural_report['try_to_skip_exons'].append(f"{exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end}|{impact}")
                            logging.info(f"add {exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end} to structural_report try_to_skip_exons")
                            break
                
                    if exon_to_skip:
                        log_key = (mrna_obj.id, mut.id, 'skip')
                        if log_key not in self.logged_splice_events:
                            logging.info(f"Simulating exon skip for {mut.id} affecting exon {exon_to_skip.start}-{exon_to_skip.end}")
                            self.logged_splice_events.add(log_key)
                        active_parts = [p for p in active_parts if not (p.start >= exon_to_skip.start and p.end <= exon_to_skip.end)]
                        structural_report['skipped_exons'].append(f"{exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end}|{impact}")
                        logging.info(f"add {exon_to_skip.chrom}:{exon_to_skip.start}-{exon_to_skip.end} to structural_report skipped_exons")

                    else:
                        logging.warning(f"Could not find a clear exon to skip for SAS mutation {mut.id} at {mut.pos}. Skipping exon skip simulation.")

                if site_type == 'donor':
                    logging.info(f"Try to simulating {mut}")
                    intron_to_retain = None
                    for intron in mrna_obj.regions['INTRON']:
                        if mrna_obj.strand == '+' and intron.start - 2 <= site_start <= intron.start + 2:
                            intron_to_retain = intron
                            structural_report['try_to_retain_introns'].append(f"{intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end}|{impact}")
                            logging.info(f"add {intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end} to structural_report try_to_retain_introns")
                            break

                        elif mrna_obj.strand == '-' and intron.end - 2 <= site_start <= intron.end + 2:
                            intron_to_retain = intron
                            structural_report['try_to_retain_introns'].append(f"{intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end}|{impact}")
                            logging.info(f"add {intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end} to structural_report try_to_retain_introns")
                            break

                    if intron_to_retain:
                        log_key = (mrna_obj.id, mut.id, 'retain')
                        if log_key not in self.logged_splice_events:
                            logging.info(f"Simulating intron retention for {mut.id} affecting intron {intron_to_retain.start}-{intron_to_retain.end}")
                            self.logged_splice_events.add(log_key)
                        upstream_exon = None
                        for exon in mrna_obj.regions['EXON']:
                            if mrna_obj.strand == '+' and exon.end + 1 == intron_to_retain.start:
                                upstream_exon = exon
                                break
                            elif mrna_obj.strand == '-' and exon.start == intron_to_retain.end + 1: 
                                upstream_exon = exon
                                break

                        if upstream_exon:
                            last_cds_idx_in_exon = -1
                            for i, part in enumerate(active_parts):
                                if part.type == 'CDS' and part.start >= upstream_exon.start and part.end <= upstream_exon.end:
                                    last_cds_idx_in_exon = i
                        
                            if last_cds_idx_in_exon != -1:
                                if mrna_obj.strand == '+':
                                    insertion_idx = last_cds_idx_in_exon + 1
                                else: 
                                    insertion_idx = last_cds_idx_in_exon
                                active_parts.insert(insertion_idx, intron_to_retain)
                                structural_report['retained_introns'].append(f"{intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end}|{impact}")
                                logging.info(f"add {intron_to_retain.chrom}:{intron_to_retain.start}-{intron_to_retain.end} to structural_report retained_introns")
                            else:
                                logging.warning(f"Could not find a matching CDS part for upstream exon {upstream_exon.start}-{upstream_exon.end} to insert retained intron.")
                        else:
                            logging.warning(f"Could not find a matching upstream exon for retained intron {intron_to_retain.start}-{intron_to_retain.end}.")

        coord_map: Dict[int, int] = {}
        part_map: Dict[int, Tuple[int, int]] = {}  
        rev_coord_map: Dict[int, int] = {} 

        current_cds_pos = 0
        temp_seq_parts = []
        
        for part_idx, part in enumerate(active_parts):
            temp_seq_parts.append(part.sequence)
            for i, g_pos in enumerate(range(part.start, part.end + 1)):
                cds_pos = current_cds_pos + i
                coord_map[g_pos] = cds_pos
                part_map[g_pos] = (part.start, part.end)
                rev_coord_map[cds_pos] = g_pos 
            current_cds_pos += len(part.sequence)
        
        cds_seq_list = list("".join(temp_seq_parts))
        retained_intron_coords = defaultdict(list)
        if structural_report.get('retained_introns'):
            for intron_str in structural_report['retained_introns']:
                try:
                    clean_str = intron_str.split('|')[0]
                    chrom, coords = clean_str.split(':')
                    start, end = map(int, coords.split('-'))
                    retained_intron_coords[chrom].append((start, end))
                    logging.info(f"retained_intron_coords: {retained_intron_coords}")
                except ValueError:
                    logging.warning(f"Could not parse retained intron coordinates: {intron_str}")
                    continue
        
        mutations_to_apply = []
        if not retained_intron_coords:
            mutations_to_apply = [m for m in mutations if 'CDS' in m.labels]
        else:
            logging.info("Retained introns found. Filtering for CDS and retained intron mutations.")
            for m in mutations:
                m_end = m.pos + len(m.ref) - 1
                is_in_retained_intron = False
                if m.chrom in retained_intron_coords:
                    for start, end in retained_intron_coords[m.chrom]:
                        if max(m.pos, start) <= min(m_end, end):
                            is_in_retained_intron = True
                            break            
                if 'CDS' in m.labels or ('INTRON' in m.labels and is_in_retained_intron):
                    mutations_to_apply.append(m)

        if not mutations_to_apply: 
            return "".join(cds_seq_list), structural_report
        sorted_muts = sorted(mutations_to_apply, key=lambda m: m.pos)
        mutation_groups = []
        if sorted_muts:
            current_group = [sorted_muts[0]]
            group_max_end = sorted_muts[0].pos + len(sorted_muts[0].ref) - 1
            for i in range(1, len(sorted_muts)):
                current_mut = sorted_muts[i]
                if current_mut.pos <= group_max_end:
                    current_group.append(current_mut)
                    group_max_end = max(group_max_end, current_mut.pos + len(current_mut.ref) - 1)
                else:
                    mutation_groups.append(current_group)
                    current_group = [current_mut]
                    group_max_end = current_mut.pos + len(current_mut.ref) - 1
            mutation_groups.append(current_group)

        def get_group_start_pos(group): return min(m.pos for m in group)

        for group in reversed(sorted(mutation_groups, key=get_group_start_pos)):
            group_mutations_to_apply = []
            sorted_muts_in_group = sorted(group, key=lambda m: m.pos)
            children_of_parent = defaultdict(list)
            is_child_variant = set()
            for i in range(len(sorted_muts_in_group)):
                for j in range(i + 1, len(sorted_muts_in_group)):
                    parent = sorted_muts_in_group[i]; child = sorted_muts_in_group[j]
                    parent_end_pos = parent.pos + len(parent.ref); child_end_pos = child.pos + len(child.ref)
                    if parent.pos <= child.pos and parent_end_pos >= child_end_pos and parent != child:
                        children_of_parent[parent.pos].append(child); is_child_variant.add(child.pos)
            for mut in sorted_muts_in_group:
                if mut.pos in is_child_variant: continue
                ref_for_replacement = mut.ref; alt_for_replacement = mut.alt
                internal_shift = 0
                if mut.pos in children_of_parent:
                    template_sequence = list(alt_for_replacement)
                    sorted_child_list = sorted(children_of_parent[mut.pos], key=lambda m: m.pos)
                    for child in sorted_child_list:
                        parent_cds_pos = coord_map.get(mut.pos); child_cds_pos = coord_map.get(child.pos)
                        if parent_cds_pos is None or child_cds_pos is None: continue
                        relative_cds_offset = child_cds_pos - parent_cds_pos
                        dynamic_internal_offset = relative_cds_offset + internal_shift
                        if dynamic_internal_offset < 0 or dynamic_internal_offset >= len(template_sequence): 
                            logging.warning(f"Child {child.id} cannot be applied; parent {mut.id} alt is too short.")
                            continue
                        current_base_on_template = "".join(template_sequence[dynamic_internal_offset : dynamic_internal_offset + len(child.ref)])
                        if current_base_on_template != child.ref:
                            logging.warning(f"Conflict: Child {child.id} expects ref '{child.ref}' but parent {mut.id} has '{current_base_on_template}' at relative pos {relative_cds_offset}.")
                            continue
                        template_sequence[dynamic_internal_offset : dynamic_internal_offset + len(child.ref)] = list(child.alt)
                        internal_shift += len(child.alt) - len(child.ref)
                    alt_for_replacement = "".join(template_sequence)
                group_mutations_to_apply.append((mut.pos, ref_for_replacement, alt_for_replacement))
            
            group_internal_shift = 0
            for pos, ref, alt in sorted(group_mutations_to_apply, key=lambda x: x[0]):
                cds_indices_to_remove = []
                first_cds_g_pos = -1

                has_cds_part = False
                has_non_cds_part = False
                ref_crosses_boundary = False
  
                affected_cds_part = []
                affected_part_indices = set()
                for i in range(len(ref)):
                    g_pos = pos + i
                    if g_pos in part_map:      
                        has_cds_part = True
                        p_start, p_end = part_map[g_pos]

                        if p_start not in affected_part_indices:
                            affected_cds_part.append((p_start, p_end))
                            affected_part_indices.add(p_start)
                    else:
                        has_non_cds_part = True
                if has_cds_part and has_non_cds_part:
                    ref_crosses_boundary = True

                affected_cds_part.sort()
                logging.info(f"pos: {pos}, ref_crosses_boundary: {ref_crosses_boundary}")

                for i in range(len(ref)):
                    g_pos = pos + i
                    if g_pos in coord_map:
                        if first_cds_g_pos == -1: first_cds_g_pos = g_pos
                        cds_indices_to_remove.append(coord_map[g_pos])
                    alt_for_cds = ""
                if ref_crosses_boundary:
                    temp_alt_slices = []
                    for p_start, p_end in affected_cds_part:
                        slice_for_this_part = self.get_alt_slice_for_cds(ref, alt, pos, p_start, p_end)
                        temp_alt_slices.append(slice_for_this_part)
                        logging.info(f"  - For part {p_start}-{p_end}, calculated slice: '{slice_for_this_part}'")
                    alt_for_cds = "".join(temp_alt_slices)
                    logging.info(f"Final combined ALT slice for CDS: '{alt_for_cds}'")
                else:
                    is_in_cds = False
                    for i in range(len(ref)):
                        if (pos + i) in coord_map:
                            is_in_cds = True
                            break
                    if is_in_cds:
                        alt_for_cds = alt 

                if not cds_indices_to_remove and not alt_for_cds:
                    continue

                num_bases_to_remove = len(cds_indices_to_remove)
                insertion_cds_pos = min(cds_indices_to_remove) if cds_indices_to_remove else coord_map.get(first_cds_g_pos)
            
                if insertion_cds_pos is not None:
                    dynamic_cds_pos = insertion_cds_pos + group_internal_shift
                    logging.info(f"cds_seq_list[dynamic_cds_pos : dynamic_cds_pos + num_bases_to_remove]: {cds_seq_list[dynamic_cds_pos : dynamic_cds_pos + num_bases_to_remove]}") 
                
                    cds_seq_list[dynamic_cds_pos : dynamic_cds_pos + num_bases_to_remove] = list(alt_for_cds) # Apply the change
                    shift = len(alt_for_cds) - num_bases_to_remove
                    group_internal_shift += shift
                    logging.info(f"alt_for_cds: {alt_for_cds}")

                logging.info(f"pos: {pos}, len(cds_seq_list_after_alt): {len(cds_seq_list)}")

        final_cds_sequence = "".join(cds_seq_list)
        return final_cds_sequence, structural_report

