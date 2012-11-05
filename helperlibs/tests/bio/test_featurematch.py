try:
    import unittest2
except ImportError:
    import unittest as unittest2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from helperlibs.bio.featurematch import FeatureMatch, find_features

class TestFeatureMatch(unittest2.TestCase):
    def setUp(self):
        self.seq = Seq("CCCAAAATGTACTCCACTATCTGCTGATTTGGG", generic_dna)
        self.feature = SeqFeature(FeatureLocation(6, 27), type="gene", strand=1)
        self.feature_seq = self.seq[3:-3]
        self.match = FeatureMatch(self.feature, self.feature_seq, 1, 3)


    def test__init(self):
        "Test FeatureMatch object creation"
        # forward strand
        m = self.match
        self.assertEqual(m.direction, "forward")
        self.assertEqual(str(m.dna), "ATGTACTCCACTATCTGCTGA")
        self.assertEqual(str(m.long_dna), str(self.feature_seq))
        self.assertEqual(str(m.promotor_region), "AAA")
        self.assertEqual(str(m.mrna), str(self.feature_seq.transcribe()))
        self.assertEqual(str(m.aas), str(self.feature.extract(self.seq).translate(to_stop=True)))

        # reverse strand
        inv_seq = self.seq.reverse_complement()
        feature_seq = inv_seq[3:-3]
        inv_feature = SeqFeature(FeatureLocation(6, 27), type="gene", strand=-1)
        m = FeatureMatch(inv_feature, feature_seq, -1, 3)
        self.assertEqual(m.direction, "reverse")
        self.assertEqual(str(m.dna), str(feature_seq[3:-3].reverse_complement()))
        self.assertEqual(str(m.long_dna), str(self.feature_seq))
        self.assertEqual(str(m.promotor_region), "AAA")
        self.assertEqual(str(m.mrna), str(feature_seq.transcribe().reverse_complement()))
        self.assertEqual(str(m.aas), str(inv_feature.extract(inv_seq).translate(to_stop=True)))


    def test_get_fasta_header(self):
        "Test FeatureMatch FASTA header creation"
        expected = ">untagged"
        self.assertEqual(self.match.get_fasta_header(), expected)

        self.match.feature.qualifiers['gene'] = ['fake']
        expected = ">fake"
        self.assertEqual(self.match.get_fasta_header(), expected)

        self.match.feature.qualifiers['locus_tag'] = ['FAKE_0001']
        expected = ">FAKE_0001"
        self.assertEqual(self.match.get_fasta_header(), expected)

        self.match.feature.qualifiers['product'] = ['Mup1']
        expected = ">FAKE_0001|Mup1"
        self.assertEqual(self.match.get_fasta_header(), expected)

        self.match.feature.qualifiers['protein_id'] = ['MUP_0001']
        expected = ">FAKE_0001|Mup1|MUP_0001"
        self.assertEqual(self.match.get_fasta_header(), expected)


    def test__str_(self):
        "Test FeatureMatch string representation"
        expected = """Feature:
	untagged
	Strand: forward
	DNA: ATGTACTCCACTATCTGCTGA
	mRNA: AAAAUGUACUCCACUAUCUGCUGAUUU
	Protein: MYSTIC"""
        self.assertMultiLineEqual(str(self.match), expected)

        self.match.feature.qualifiers['gene'] = ["fake"]
        expected = """Feature:
	Tag: fake
	Strand: forward
	DNA: ATGTACTCCACTATCTGCTGA
	mRNA: AAAAUGUACUCCACUAUCUGCUGAUUU
	Protein: MYSTIC"""
        self.assertMultiLineEqual(str(self.match), expected)

        self.match.feature.qualifiers['locus_tag'] = ['FAKE_0001']
        expected = """Feature:
	Tag: FAKE_0001
	Strand: forward
	DNA: ATGTACTCCACTATCTGCTGA
	mRNA: AAAAUGUACUCCACUAUCUGCUGAUUU
	Protein: MYSTIC"""
        self.assertMultiLineEqual(str(self.match), expected)

        self.match.feature.qualifiers['product'] = ['Mup1']
        expected = """Feature:
	Tag: FAKE_0001
	Strand: forward
	Product: Mup1
	DNA: ATGTACTCCACTATCTGCTGA
	mRNA: AAAAUGUACUCCACUAUCUGCUGAUUU
	Protein: MYSTIC"""
        self.assertMultiLineEqual(str(self.match), expected)

        self.match.feature.qualifiers['protein_id'] = ['MUP_0001']
        expected = """Feature:
	Tag: FAKE_0001
	Strand: forward
	Product: Mup1
	Protein ID: MUP_0001
	DNA: ATGTACTCCACTATCTGCTGA
	mRNA: AAAAUGUACUCCACUAUCUGCUGAUUU
	Protein: MYSTIC"""
        self.assertMultiLineEqual(str(self.match), expected)


    def test_dna_fasta(self):
        "Test FeatureMatch DNA FASTA output"
        expected = "%s\n%s" % (self.match.get_fasta_header(), self.match.dna)
        self.assertMultiLineEqual(self.match.dna_fasta(), expected)


    def test_long_dna_fasta(self):
        "Test FeatureMatch long DNA FASTA output"
        expected = "%s\n%s" % (self.match.get_fasta_header(), self.match.long_dna)
        self.assertMultiLineEqual(self.match.long_dna_fasta(), expected)


    def test_mrna_fasta(self):
        "Test FeatureMatch mRNA FASTA output"
        expected = "%s\n%s" % (self.match.get_fasta_header(), self.match.mrna)
        self.assertMultiLineEqual(self.match.mrna_fasta(), expected)


    def test_protein_fasta(self):
        "Test FeatureMatch protein FASTA output"
        expected = "%s\n%s" % (self.match.get_fasta_header(), self.match.aas)
        self.assertMultiLineEqual(self.match.protein_fasta(), expected)


    def test_promotor_fasta(self):
        "Test FeatureMatch promotor DNA FASTA output"
        expected = "%s\n%s" % (self.match.get_fasta_header(), self.match.promotor_region)
        self.assertMultiLineEqual(self.match.promotor_fasta(), expected)


class TestFindFeatures(unittest2.TestCase):
    def setUp(self):
        seq = Seq("CCCAAAATGTACTCCACTATCTGCTGATTTGGG", generic_dna)
        feature = SeqFeature(FeatureLocation(6, 27), type="CDS", strand=1)
        self.record = SeqRecord(seq, features=[feature])


    def test_find_features(self):
        "Test find_features convenience function"
        found = find_features([self.record], utr_len=3)
        self.assertEqual(len(found), 1)
        self.assertMultiLineEqual(str(found[0].feature), str(self.record.features[0]))
