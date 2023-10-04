"""Kakapo workflow: Find ORFs."""

import json
from functools import partial
from os.path import join as opj
from sys import exit

from kakapo.data.start_codon_context import contexts as atg_contexts
from kakapo.tools.bioio import (seq_records_to_dict,
                                trim_desc_to_first_space_in_fasta_text,
                                write_fasta)
from kakapo.tools.blast import collate_blast_results
from kakapo.tools.orf import find_orf_for_blast_hit
from kakapo.tools.seq import SEQ_TYPE_DNA, reverse_complement, translate
from kakapo.tools.seqtk import seqtk_extract_reads
from kakapo.utils.logging import Log


def find_orfs_translate(ss, assemblies, dir_prj_transcripts, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod, tax, tax_group, tax_ids_user,
                        min_overlap, ss_organelle):

    if len(assemblies) > 0:
        if seqtk is None:
            Log.err('seqtk is not available. Cannot continue. Exiting.', '')
            exit(0)

    for a in assemblies:

        if ('blast_hits_aa__' + ss) not in a:
            continue

        assmbl_name = a['name']
        tax_id = a['tax_id']

        parsed_hits = a['blast_hits_aa__' + ss]

        a_path = a['path']

        gc_tt = a['gc_tt']
        if tax.is_eukaryote(tax_id) is True:
            if ss_organelle == 'mitochondrion':
                gc_tt = a['gc_tt_mito']
            if tax.contains_plastid(tax_id) is True:
                if ss_organelle == 'plastid':
                    gc_tt = a['gc_tt_plastid']

        transcripts_nt_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt__' + ss
            + '.fasta')

        transcripts_nt_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt_orf__' + ss
            + '.fasta')

        transcripts_aa_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_aa_orf__' + ss
            + '.fasta')

        transcripts_nt = {}
        transcripts_nt_orf = {}
        transcripts_aa_orf = {}

        transcripts_with_acceptable_orfs = set()

        ann_key = 'annotations__'

        a[ann_key + ss] = {}

        collated = collate_blast_results(parsed_hits)

        ######################################################################
        # Use seqtk to sample the assembly FASTA file for sequences with
        # BLAST hits. This increases the speed substantially when the assembly
        # file is large.
        temp_a_file = opj(dir_temp, 'temp__' + ss + '.fasta')
        temp_s_file = opj(dir_temp, 'temp__' + ss + '.txt')
        sseqids_subsample = []
        for hit in collated:
            target_name = hit['sseqid']
            sseqids_subsample.append(target_name)
        sseqids_subsample_text = '\n'.join(sseqids_subsample)
        with open(temp_s_file, 'w') as f:
            f.write(sseqids_subsample_text)
        seqtk_extract_reads(seqtk,
                            in_file=a_path,
                            out_file=temp_a_file,
                            ids_file=temp_s_file)

        with open(temp_a_file, 'r') as f:
            _ = f.read()

        if _.strip() == '':
            continue

        print()
        Log.log(inf='Analyzing BLAST hits', s='=' * 113 + '\n')
        Log.log(msg='Assembly:', s=assmbl_name, timestamp=False)
        Log.log(msg='Search Strategy:', s=ss + '\n\n' + '-' * 134 + '\n', timestamp=False)

        parsed_fasta = trim_desc_to_first_space_in_fasta_text(_, SEQ_TYPE_DNA)
        parsed_fasta = seq_records_to_dict(parsed_fasta)
        ######################################################################

        all_kakapo_results = {}
        json_dump_file_path = opj(dir_prj_transcripts, assmbl_name
                                  + '_ann_kakapo__' + ss + '.json')

        for hit in collated:

            target_name = hit['sseqid']
            target_seq = parsed_fasta[target_name]
            query_name = hit['qseqid']
            hit_evalue = hit['evalue']

            # Prepend assembly name to the sequence name:
            if prepend_assmbl is True:
                target_name = assmbl_name + '__' + target_name
                # Also prepend taxonomic info to the sequence name:
                if tax_id is not None:
                    fm = tax.higher_rank_for_taxid(tax_id, rank='family')
                    if fm is not None:
                        target_name = fm + '__' + target_name

            hit_start = hit['start']
            hit_end = hit['end']
            hit_frame = hit['frame']

            if allow_non_aug is True:
                start_codons = gc_tt.start_codons_ambiguous
            else:
                start_codons = ['ATG']

            stop_codons = gc_tt.stop_codons_ambiguous

            ##################################################################
            if tax_id is not None:
                tax_ids_for_orf = (tax_id, )
            else:
                tax_ids_for_orf = tax_ids_user

            cntx_txids_avail = tuple(
                sorted(set(map(lambda x: int(x.split('_')[0]),
                               atg_contexts.keys()))))

            cntx_taxid = set()
            for txid in tax_ids_for_orf:
                tax_path = partial(tax.path_between_taxids, txid)
                path_len = tuple(map(len,
                                     tuple(map(tax_path, cntx_txids_avail))))
                cntx_taxid.add(cntx_txids_avail[path_len.index(min(path_len))])
            cntx_taxid = tuple(cntx_taxid)[0]

            cntx_l_key = str(cntx_taxid) + '_L'
            cntx_r_key = str(cntx_taxid) + '_R'

            cntx_l = atg_contexts[cntx_l_key]
            cntx_r = atg_contexts[cntx_r_key]
            ##################################################################

            orf_log_str = ('grade'.rjust(5) + 'ovrlp'.rjust(7)
                           + 'cntx'.rjust(6) + 'length'.center(9)
                           + 'cntx_l'.rjust(8) + 'cntx_r'.rjust(15) + '\n')

            orf = find_orf_for_blast_hit(
                seq=target_seq,
                frame=hit_frame,
                hit_start=hit_start,
                hit_end=hit_end,
                stop_codons=stop_codons,
                start_codons=start_codons,
                context_l=cntx_l,
                context_r=cntx_r,
                min_overlap=min_overlap,
                min_len=min_target_orf_len,
                max_len=max_target_orf_len,
                allow_no_strt_cod=allow_no_strt_cod,
                allow_no_stop_cod=allow_no_stop_cod)

            orf_log_str += orf[2]

            rev_comp_def_str = ''
            if hit_frame > 0:
                ann_hit_b = hit_start
                ann_hit_e = hit_end
            else:
                target_seq = reverse_complement(target_seq)
                ann_hit_b = len(target_seq) - hit_start
                ann_hit_e = len(target_seq) - hit_end
                rev_comp_def_str = '; RevComp'

            target_def = target_name + ' ' + query_name + rev_comp_def_str

            a[ann_key + ss][target_name] = {}

            good_orfs = orf[0]
            bad_orfs = orf[1]

            if len(good_orfs) > 0:
                a[ann_key + ss][target_name]['orfs_good'] = dict()
                orfs_good_dict = a[ann_key + ss][target_name]['orfs_good']
                orf_log_str += '\n' + 'VALID ' + '-' * 128 + '\n'

                for i, good_orf in enumerate(good_orfs):

                    good_orf_frame = good_orf[2]

                    if good_orf_frame > 0:
                        ann_orf_b = good_orf[0]
                        ann_orf_e = good_orf[1] + 3
                        orf_seq = target_seq[ann_orf_b:ann_orf_e]
                    else:
                        ann_orf_b = len(target_seq) - good_orf[1]
                        ann_orf_e = len(target_seq) - good_orf[0] + 3
                        orf_seq = target_seq[ann_orf_b:ann_orf_e]

                    orf_good_dict = dict()
                    orf_good_dict['orf_begin'] = ann_orf_b
                    orf_good_dict['orf_end'] = ann_orf_e
                    orf_good_dict['orf_frame'] = abs(good_orf_frame)
                    orf_good_dict['orf_grade'] = good_orf[3]
                    orf_good_dict['orf_tt_id'] = str(gc_tt.gc_id)
                    orf_good_dict['orf_tt_name'] = gc_tt.gc_name

                    orfs_good_dict['ORF{:03d}'.format(i + 1)] = orf_good_dict

                    target_def_orf = (target_name
                                      + '__ORF{:03d}'.format(i + 1) + ' '
                                      + query_name + rev_comp_def_str)

                    transcripts_nt_orf[target_def_orf] = orf_seq

                    transcripts_with_acceptable_orfs.add(target_name)

                    transl_seq = translate(orf_seq,
                                           gc_tt.table_ambiguous,
                                           start_codons)

                    transcripts_aa_orf[target_def_orf] = transl_seq[:-1]

            else:
                orf_log_str += '\n' + 'NOT VALID ' + '-' * 124 + '\n'

            Log.log(msg='Transcript:', s=target_name, timestamp=False)
            Log.log(msg='     Query:', s=query_name + '\n\n' + orf_log_str, timestamp=False)

            if len(bad_orfs) > 0:
                a[ann_key + ss][target_name]['orfs_bad'] = dict()
                orfs_bad_dict = a[ann_key + ss][target_name]['orfs_bad']

                for i, bad_orf in enumerate(bad_orfs):

                    bad_orf_frame = bad_orf[2]

                    if bad_orf_frame > 0:
                        ann_orf_b = bad_orf[0]
                        ann_orf_e = bad_orf[1] + 3
                        orf_seq = target_seq[ann_orf_b:ann_orf_e]
                    else:
                        ann_orf_b = len(target_seq) - bad_orf[1]
                        ann_orf_e = len(target_seq) - bad_orf[0] + 3
                        orf_seq = target_seq[ann_orf_b:ann_orf_e]

                    orf_bad_dict = dict()
                    orf_bad_dict['orf_begin'] = ann_orf_b
                    orf_bad_dict['orf_end'] = ann_orf_e
                    orf_bad_dict['orf_frame'] = abs(bad_orf_frame)
                    orf_bad_dict['orf_grade'] = bad_orf[3]
                    orf_bad_dict['orf_tt_id'] = str(gc_tt.gc_id)
                    orf_bad_dict['orf_tt_name'] = gc_tt.gc_name

                    orfs_bad_dict['ORF{:03d}'.format(i + 1)] = orf_bad_dict

            transcripts_nt[target_def] = target_seq

            a[ann_key + ss][target_name]['blast_hit'] = dict()
            blast_hit_dict = a[ann_key + ss][target_name]['blast_hit']
            blast_hit_dict['query_name'] = query_name
            blast_hit_dict['query_id'] = ss
            blast_hit_dict['evalue'] = hit_evalue
            blast_hit_dict['frame'] = abs(hit_frame)
            blast_hit_dict['blast_hit_begin'] = ann_hit_b
            blast_hit_dict['blast_hit_end'] = ann_hit_e

            # Collect ORF and BLAST hit annotations for downstream use. ######
            kakapo_json = [{}]
            kakapo_json[0]['kakapo_annotations__' + ss] = (
                a[ann_key + ss][target_name])
            all_kakapo_results[target_name] = kakapo_json
            ##################################################################

        # --------------------------------------------------------------------

        Log.log(msg='Assembly:', s=assmbl_name, timestamp=False)
        Log.log(msg='Search Strategy:', s=ss, timestamp=False)
        Log.log(msg='Transcripts:', s=str(len(transcripts_nt)), timestamp=False)
        Log.log(msg='Transcripts with acceptable ORFs:',
                s=str(len(transcripts_with_acceptable_orfs)) + '\n'
                + '=' * 134, timestamp=False)

        if len(transcripts_nt) > 0:
            write_fasta(transcripts_nt, transcripts_nt_fasta_file)
            a['transcripts_nt_fasta_file__' + ss] = transcripts_nt_fasta_file
        else:
            a['transcripts_nt_fasta_file__' + ss] = None

        if len(transcripts_nt_orf) > 0:
            write_fasta(transcripts_nt_orf, transcripts_nt_orf_fasta_file)
            a['transcripts_nt_orf_fasta_file__' + ss] = transcripts_nt_orf_fasta_file
        else:
            a['transcripts_nt_orf_fasta_file__' + ss] = None

        if len(transcripts_aa_orf) > 0:
            write_fasta(transcripts_aa_orf, transcripts_aa_orf_fasta_file)
            a['transcripts_aa_orf_fasta_file__' + ss] = transcripts_aa_orf_fasta_file
        else:
            a['transcripts_aa_orf_fasta_file__' + ss] = None

        # Save ORF and BLAST hit annotations for downstream use.--------------
        with open(json_dump_file_path, 'w') as f:
            json.dump(all_kakapo_results, f, sort_keys=True, indent=4)
        # --------------------------------------------------------------------
