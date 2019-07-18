# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from os.path import join as opj

import unittest
from datetime import datetime

from kakapo.seq import Seq
from kakapo.seq import SeqRecord

from kakapo.seq import NTSeq
from kakapo.seq import DNASeq
from kakapo.seq import RNASeq
from kakapo.seq import AASeq
from kakapo.seq import reverse_complement

from kakapo.seq import SEQ_TYPE_NT
from kakapo.seq import SEQ_TYPE_DNA
from kakapo.seq import SEQ_TYPE_RNA
from kakapo.seq import SEQ_TYPE_AA

from kakapo.seq import SEQ_TYPES
from kakapo.seq import MOL_TO_SEQ_TYPE_MAP

from kakapo.bioio import _parse_gbseq_xml_text

from tests import test_data_dir_path


def setUpModule():
    print('\nsetUpModule kakapoSeqTests')


def tearDownModule():
    print('\n\ntearDownModule kakapoSeqTests')


class kakapoSeqTests(unittest.TestCase):

    print('kakapoSeqTests')

    def test_reverse_complement(self):

        print('\ntest_reverse_complement')

        seq = 'ACGTRYMKWSBDHV'
        rev_com_seq = 'BDHVSWMKRYACGT'
        test_result = reverse_complement(seq)
        self.assertEqual(test_result, rev_com_seq)

    def test_Seq(self):

        print('\ntest_Seq')

        seq = Seq(seq='ACGTUBDHKMNRSVWY', seq_type=SEQ_TYPE_NT)
        self.assertTrue(isinstance(seq, NTSeq))
        self.assertEqual(seq.length, 16)

        seq = Seq(seq='ACGTBDHKMNRSVWY', seq_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq, DNASeq))
        self.assertEqual(seq.length, 15)

        seq = Seq(seq='ACGUBDHKMNRSVWY', seq_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq, RNASeq))
        self.assertEqual(seq.length, 15)

        seq = Seq(seq='GAVLIPFYCMHKRWSTDENQ.BZX', seq_type=SEQ_TYPE_AA)
        self.assertTrue(isinstance(seq, AASeq))
        self.assertEqual(seq.length, 24)

    def test_SeqRecord(self):

        print('\ntest_SeqRecord')

        seq_record = SeqRecord(seq='ACGTUBDHKMNRSVWY', mol_type=SEQ_TYPE_NT)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, NTSeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWY', mol_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWYU', mol_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWYT', mol_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='GAVLIPFYCMHKRWSTDENQ.BZX',
                               mol_type=SEQ_TYPE_AA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, AASeq))
        self.assertEqual(seq_record.seq.length, 24)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWY', mol_type='DNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='mRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='rRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='tRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='GAVLIPFYCMHKRWSTDENQ.BZX', mol_type='AA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, AASeq))
        self.assertEqual(seq_record.seq.length, 24)

        with open(opj(test_data_dir_path, 'gbseq_xml_sample.xml'), 'r') as f:
            results = f.read()

        results = _parse_gbseq_xml_text(results)

        self.assertEqual(len(results), 6)

        for r in results:

            self.assertEqual(len(r.keys()), 14)

            seq_record = SeqRecord(
                seq=r['seq'],
                mol_type=r['mol_type'],
                accession=r['accession'],
                version=r['version'],
                description=r['definition'],
                strandedness=r['strandedness'],
                topology=r['topology'],
                division=r['division'],
                date_create=(datetime.strptime(
                    r['date_create'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
                date_update=(datetime.strptime(
                    r['date_update'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
                taxid=r['taxid'],
                organism=r['organism'],
                features=r['features']
            )

            self.assertTrue(isinstance(seq_record, SeqRecord))

            seq_type = None

            if r['mol_type'] in SEQ_TYPES:
                seq_type = r['mol_type']
            elif r['mol_type'] in MOL_TO_SEQ_TYPE_MAP:
                seq_type = MOL_TO_SEQ_TYPE_MAP[r['mol_type']]

            if seq_type is SEQ_TYPE_NT:
                self.assertTrue(isinstance(seq_record.seq, NTSeq))
            elif seq_type is SEQ_TYPE_DNA:
                self.assertTrue(isinstance(seq_record.seq, DNASeq))
            elif seq_type is SEQ_TYPE_RNA:
                self.assertTrue(isinstance(seq_record.seq, RNASeq))
            elif seq_type is SEQ_TYPE_AA:
                self.assertTrue(isinstance(seq_record.seq, AASeq))
