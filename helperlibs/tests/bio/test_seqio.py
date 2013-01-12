try:
    import unittest2
except ImportError:
    import unittest as unittest2
import Bio.SeqIO
from helperlibs.bio import seqio
from minimock import TraceTracker, assert_same_trace, mock, restore


class DummyHandle(object):
    def __init__(self, name):
        self.name = name


    def __repr__(self):
        return "DummyHandle(%r)" % self.name


class TestSeqIO(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        self.handle = DummyHandle("test.gbk")


    def tearDown(self):
        restore()


    def test__get_seqtype_from_ext(self):
        "Test guessing the sequence type from the file extension for handle input"
        gbk_h = DummyHandle("test.gbk")
        gb_h = DummyHandle("test.gb")
        genbank_h = DummyHandle("test.genbank")
        gbff_h = DummyHandle("test.gbff")
        embl_h = DummyHandle("test.embl")
        emb_h = DummyHandle("test.emb")
        fa_h = DummyHandle("test.fa")
        fasta_h = DummyHandle("test.fasta")
        fna_h = DummyHandle("test.fna")
        faa_h = DummyHandle("test.faa")
        fas_h = DummyHandle("test.fas")
        invalid_h = DummyHandle("test.invalid")

        for handle in (gbk_h, gb_h, genbank_h, gbff_h):
            self.assertEqual("genbank", seqio._get_seqtype_from_ext(handle))

        for handle in (embl_h, emb_h):
            self.assertEqual("embl", seqio._get_seqtype_from_ext(handle))

        for handle in (fa_h, fasta_h, fna_h, faa_h, fas_h):
            self.assertEqual("fasta", seqio._get_seqtype_from_ext(handle))

        self.assertRaises(ValueError, seqio._get_seqtype_from_ext, invalid_h)


    def test__get_seqtype_from_ext_string(self):
        "Test guessing the sequence type from the file extension for string input"
        for string in ("test.gbk", "test.gb", "test.genbank", "test.gbff"):
            self.assertEqual("genbank", seqio._get_seqtype_from_ext(string))

        for string in ("test.embl", "test.emb"):
            self.assertEqual("embl", seqio._get_seqtype_from_ext(string))

        for string in ("test.fa", "test.fasta", "test.fna", "test.faa", "test.fas"):
            self.assertEqual("fasta", seqio._get_seqtype_from_ext(string))

        self.assertRaises(ValueError, seqio._get_seqtype_from_ext, "test.invalid")


    def test_parse(self):
        "Test running the Bio.SeqIO parser"
        mock("Bio.SeqIO.parse", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.parse(DummyHandle('test.gbk'), 'genbank')"
        seqio.parse(self.handle)
        assert_same_trace(self.tt, expected_trace)


    def test_read(self):
        "Test reading a single sequence via Bio.SeqIO"
        mock("Bio.SeqIO.read", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.read(DummyHandle('test.gbk'), 'genbank')"
        seqio.read(self.handle)
        assert_same_trace(self.tt, expected_trace)


    def test_write(self):
        "Test writing Bio.SeqIO records"
        mock("Bio.SeqIO.write", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.write(['fake'], DummyHandle('test.gbk'), 'genbank')"
        seqio.write(['fake'], self.handle, "genbank")
        assert_same_trace(self.tt, expected_trace)
