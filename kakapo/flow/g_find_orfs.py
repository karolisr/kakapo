# -*- coding: utf-8 -*-

"""Kakapo workflow: Find ORFs."""

import json

from functools import partial
from io import StringIO
from os.path import join as opj
from sys import exit

from kakapo.bioio import read_fasta
from kakapo.bioio import trim_desc_to_first_space_in_fasta_text
from kakapo.bioio import write_fasta
from kakapo.blast import collate_blast_results
from kakapo.data.start_codon_context import contexts as atg_contexts
from kakapo.orf import find_orf_for_blast_hit
from kakapo.seq import reverse_complement, translate
from kakapo.seqtk import seqtk_extract_reads


def find_orfs_translate(ss, assemblies, dir_prj_transcripts, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod, tax, tax_group, tax_ids_user,
                        min_overlap, organelle, linfo=print):

    if len(assemblies) > 0:
        linfo('Analyzing BLAST hits for assemblies [' + ss + ']')
        if seqtk is None:
            linfo('seqtk is not available. ' +
                  'Cannot continue. Exiting.')
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
            if organelle == 'mitochondrion':
                gc_tt = a['gc_tt_mito']
            if tax.contains_plastid(tax_id) is True:
                if organelle == 'plastid':
                    gc_tt = a['gc_tt_plastid']

        transcripts_nt_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt__' + ss + '.fasta')

        transcripts_nt_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt_orf__' + ss + '.fasta')

        transcripts_aa_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_aa_orf__' + ss + '.fasta')

        transcripts_nt = {}
        transcripts_nt_orf = {}
        transcripts_aa_orf = {}

        a['annotations__' + ss] = {}

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

        linfo(assmbl_name)

        _ = trim_desc_to_first_space_in_fasta_text(_)

        parsed_fasta = read_fasta(StringIO(_))
        ######################################################################

        all_kakapo_results = {}
        json_dump_file_path = opj(dir_prj_transcripts, assmbl_name +
                                  '_ann_kakapo__' + ss + '.json')

        if len(collated) > 0:
            print('\n' + '-' * 90 + '\n')

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

            orf_log_str = target_name.center(90) + '\n'
            orf_log_str += query_name.center(90) + '\n\n'

            orf_log_str += ('grade'.rjust(6) + 'ovrlp'.rjust(7) +
                            'cntx'.rjust(6) + 'length'.center(9) +
                            'cntx_l'.rjust(7) + 'cntx_r'.rjust(15) + '\n')

            orf = find_orf_for_blast_hit(
                seq=target_seq,
                frame=hit_frame,
                hit_start=hit_start,
                hit_end=hit_end,
                stop_codons=stop_codons,
                start_codons=start_codons,
                context_l=cntx_l,
                context_r=cntx_r,
                min_overlap=min_overlap)

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

            a['annotations__' + ss][target_name] = {}

            good_orf = orf[0]
            bad_orfs = orf[1]

            if good_orf is not None:

                if hit_frame > 0:
                    ann_orf_b = good_orf[0]
                    ann_orf_e = good_orf[1] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]
                else:
                    ann_orf_b = len(target_seq) - good_orf[1]
                    ann_orf_e = len(target_seq) - good_orf[0] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]

                ##############################################################
                valid_orf = True
                invalid_orf_reason = ''

                if allow_non_aug is False and \
                        orf_seq[0:3] != 'ATG':
                    valid_orf = False
                    invalid_orf_reason = 'Start codon is not ATG.'

                elif allow_no_strt_cod is False and \
                        orf_seq[0:3] not in start_codons:
                    valid_orf = False
                    invalid_orf_reason = 'No start codon.'

                elif allow_no_stop_cod is False and \
                        orf_seq[-3:] not in stop_codons:
                    valid_orf = False
                    invalid_orf_reason = 'No stop codon.'

                elif len(orf_seq) < min_target_orf_len:
                    valid_orf = False
                    invalid_orf_reason = 'ORF is not long enough.'

                elif len(orf_seq) > max_target_orf_len:
                    valid_orf = False
                    invalid_orf_reason = 'ORF is too long.'
                ##############################################################

                if valid_orf is True:

                    orf_log_str += '\n' + 'VALID'.center(90) + '\n'

                    a['annotations__' + ss][target_name]['orf_begin'] = ann_orf_b
                    a['annotations__' + ss][target_name]['orf_end'] = ann_orf_e
                    a['annotations__' + ss][target_name]['orf_grade'] = good_orf[3]
                    a['annotations__' + ss][target_name]['orf_tt_id'] = str(gc_tt.gc_id)
                    a['annotations__' + ss][target_name]['orf_tt_name'] = gc_tt.gc_name

                    transcripts_nt_orf[target_def] = orf_seq

                    transl_seq = translate(orf_seq,
                                           gc_tt.table_ambiguous,
                                           start_codons)

                    transcripts_aa_orf[target_def] = transl_seq[:-1]

                else:
                    msg = 'INVALID: ' + invalid_orf_reason
                    orf_log_str += '\n' + msg.center(90) + '\n'
                    bad_orfs.append(good_orf)

            else:
                orf_log_str += '\n' + 'INVALID'.center(90) + '\n'

            orf_log_str += '\n' + '-' * 90 + '\n'
            print(orf_log_str)

            if len(bad_orfs) > 0:
                a['annotations__' + ss][target_name]['orfs_bad'] = list()
                orfs_bad_list = a['annotations__' + ss][target_name]['orfs_bad']

            for bad_orf in bad_orfs:

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

                orfs_bad_list.append(orf_bad_dict)

            transcripts_nt[target_def] = target_seq

            a['annotations__' + ss][target_name]['query_name'] = query_name
            a['annotations__' + ss][target_name]['evalue'] = hit_evalue
            a['annotations__' + ss][target_name]['frame'] = abs(hit_frame)
            a['annotations__' + ss][target_name]['blast_hit_begin'] = ann_hit_b
            a['annotations__' + ss][target_name]['blast_hit_end'] = ann_hit_e

            # Collect ORF and BLAST hit annotations for downstream use. ######
            kakapo_json = [{}]
            kakapo_json[0]['kakapo_annotations__' + ss] = (
                a['annotations__' + ss][target_name])
            all_kakapo_results[target_name] = kakapo_json
            ##################################################################

        # --------------------------------------------------------------------

        linfo('Transcripts: ' + str(len(transcripts_nt)))

        if len(transcripts_nt) > 0:
            write_fasta(transcripts_nt, transcripts_nt_fasta_file)
            a['transcripts_nt_fasta_file__' + ss] = transcripts_nt_fasta_file
        else:
            a['transcripts_nt_fasta_file__' + ss] = None

        linfo('Transcripts with acceptable ORFs: ' +
              str(len(transcripts_nt_orf)))

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
