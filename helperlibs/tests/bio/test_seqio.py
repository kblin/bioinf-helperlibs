try:
    import unittest2
except ImportError:
    import unittest as unittest2
import Bio.SeqIO
from helperlibs.bio import seqio
from helperlibs.tests import get_file_path
from minimock import TraceTracker, assert_same_trace, mock, restore


class DummyHandle(object):
    def __init__(self, name):
        self.name = name


    def __repr__(self):
        return "DummyHandle(%r)" % self.name


class TestSeqIOSeqtype(unittest2.TestCase):
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
        gbk_gz_h = DummyHandle("test.gbk.gz")
        gb_gz_h = DummyHandle("test.gb.gz")
        genbank_gz_h = DummyHandle("test.genbank.gz")
        gbff_gz_h = DummyHandle("test.gbff.gz")
        invalid_h = DummyHandle("test.invalid")

        for handle in (gbk_h, gb_h, genbank_h, gbff_h):
            self.assertEqual("genbank", seqio._get_seqtype_from_ext(handle))

        for handle in (embl_h, emb_h):
            self.assertEqual("embl", seqio._get_seqtype_from_ext(handle))

        for handle in (fa_h, fasta_h, fna_h, faa_h, fas_h):
            self.assertEqual("fasta", seqio._get_seqtype_from_ext(handle))

        for handle in (gbk_gz_h, gb_gz_h, genbank_gz_h, gbff_gz_h):
            self.assertEqual("gz-genbank", seqio._get_seqtype_from_ext(handle))

        self.assertRaises(ValueError, seqio._get_seqtype_from_ext, invalid_h)


    def test__get_seqtype_from_ext_string(self):
        "Test guessing the sequence type from the file extension for string input"
        for string in ("test.gbk", "test.gb", "test.genbank", "test.gbff"):
            self.assertEqual("genbank", seqio._get_seqtype_from_ext(string))

        for string in ("test.embl", "test.emb"):
            self.assertEqual("embl", seqio._get_seqtype_from_ext(string))

        for string in ("test.fa", "test.fasta", "test.fna", "test.faa", "test.fas"):
            self.assertEqual("fasta", seqio._get_seqtype_from_ext(string))

        for string in ("test.gbk.gz", "test.gb.gz", "test.genbank.gz", "test.gbff.gz"):
            self.assertEqual("gz-genbank", seqio._get_seqtype_from_ext(string))

        self.assertRaises(ValueError, seqio._get_seqtype_from_ext, "test.invalid")


    def test__guess_seqtype_from_file_genbank_correct(self):
        "Test guessing the sequence type from correct genbank contents"
        with open(get_file_path('melanin.gbk'), 'rU') as h:
            self.assertEqual("genbank", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("genbank", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_embl_corect(self):
        "Test guessing the sequence type from correct embl contents"
        with open(get_file_path('melanin.embl'), 'rU') as h:
            self.assertEqual("embl", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("embl", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_fasta_correct(self):
        "Test guessing the sequence type from correct fasta contents"
        with open(get_file_path('melanin.fasta'), 'rU') as h:
            self.assertEqual("fasta", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("fasta", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_genbank_no_header(self):
        "Test guessing the sequence type from a genbank file without header"
        with open(get_file_path('no_header.gbk'), 'rU') as h:
            self.assertEqual("genbank", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("genbank", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_embl_no_header(self):
        "Test guessing the sequence type from an embl file without header"
        with open(get_file_path('no_header.embl'), 'rU') as h:
            self.assertEqual("embl", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("embl", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_fasta_no_header(self):
        "Test guessing the sequence type from a fasta file without header"
        with open(get_file_path('no_header.fasta'), 'rU') as h:
            self.assertEqual("fasta", seqio._guess_seqtype_from_file(h))
            h.seek(0)
            string_seq = h.read()
            self.assertEqual("fasta", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_fasta_no_header_lower_case(self):
        "Test guessing the sequence type from a lower case fasta file without header"
        with open(get_file_path('no_header.fasta'), 'rU') as h:
            string_seq = h.read().lower()
            self.assertEqual("fasta", seqio._guess_seqtype_from_file(string_seq))


    def test__guess_seqtype_from_file_raises_error(self):
        "Test guessing the sequence type from file raises error when it fails"
        self.assertRaises(ValueError, seqio._guess_seqtype_from_file, 'bad & invalid')


class TestSeqIODummy(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        self.handle = DummyHandle("test.gbk")


    def tearDown(self):
        restore()


    def test_parse_calls_biopython(self):
        "Test running the Bio.SeqIO parser"
        mock("Bio.SeqIO.parse", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.parse(DummyHandle('test.gbk'), 'genbank')"
        seqio.parse(self.handle)
        assert_same_trace(self.tt, expected_trace)


    def test_read_calls_biopython(self):
        "Test reading a single sequence via Bio.SeqIO"
        mock("Bio.SeqIO.read", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.read(DummyHandle('test.gbk'), 'genbank')"
        seqio.read(self.handle)
        assert_same_trace(self.tt, expected_trace)


    def test_write_calls_biopython(self):
        "Test writing Bio.SeqIO records"
        mock("Bio.SeqIO.write", tracker=self.tt, returns=[])
        expected_trace = "    Called Bio.SeqIO.write(['fake'], DummyHandle('test.gbk'), 'genbank')"
        seqio.write(['fake'], self.handle, "genbank")
        assert_same_trace(self.tt, expected_trace)


class TestSeqIORobust(unittest2.TestCase):
    def test_parse_genbank_valid(self):
        "Test parsing a valid genbank record"
        with open(get_file_path('melanin.gbk'), 'rU') as h:
            records = list(seqio.parse(h))
        self.assertEqual(1, len(records))


    def test_parse_embl_valid(self):
        "Test parsing a valid embl record"
        with open(get_file_path('melanin.embl'), 'rU') as h:
            records = list(seqio.parse(h))
        self.assertEqual(1, len(records))


    def test_parse_fasta_valid(self):
        "Test parsing a valid fasta record"
        with open(get_file_path('melanin.fasta'), 'rU') as h:
            records = list(seqio.parse(h))
        self.assertEqual(1, len(records))


    def test_parse_genbank_no_header(self):
        "Test parsing a genbank record without header"
        with open(get_file_path('no_header.gbk'), 'rU') as h:
            # plain BioPython parsing should fail
            records = list(seqio.parse(h))
            self.assertEqual(0, len(records))
            h.seek(0)
            # robust parsing should work
            records = list(seqio.parse(h, robust=True))
            self.assertEqual(1, len(records))


    def test_parse_embl_no_header(self):
        "Test parsing an embl record without header"
        with open(get_file_path('no_header.embl'), 'rU') as h:
            # plain BioPython parsing should fail
            records = list(seqio.parse(h))
            self.assertEqual(0, len(records))
            h.seek(0)
            # robust parsing should work
            records = list(seqio.parse(h, robust=True))
            self.assertEqual(1, len(records))


    def test_parse_fasta_no_header(self):
        "Test parsing a fasta record without header"
        with open(get_file_path('no_header.fasta'), 'rU') as h:
            # plain BioPython parsing should fail
            records = list(seqio.parse(h))
            self.assertEqual(0, len(records))
            h.seek(0)
            # robust parsing should work
            records = list(seqio.parse(h, robust=True))
            self.assertEqual(1, len(records))


    def test_read_genbank_valid(self):
        "Test reading a valid genbank record"
        with open(get_file_path('melanin.gbk'), 'rU') as h:
            record = seqio.read(h)
        self.assertEqual("AB070938.1", record.id)


    def test_read_embl_valid(self):
        "Test reading a valid embl record"
        with open(get_file_path('melanin.embl'), 'rU') as h:
            record = seqio.read(h)
        self.assertEqual("AB070938.1", record.id)


    def test_read_fasta_valid(self):
        "Test reading a valid fasta record"
        with open(get_file_path('melanin.fasta'), 'rU') as h:
            record = seqio.read(h)
        self.assertEqual("AB070938", record.id)


    def test_read_genbank_no_header(self):
        "Test reading a genbank record without header"
        with open(get_file_path('no_header.gbk'), 'rU') as h:
            # plain BioPython reading should fail
            self.assertRaises(ValueError, seqio.read, h)
            h.seek(0)
            # robust reading should work
            record = seqio.read(h, robust=True)
            self.assertEqual("DUMMY", record.id)


    def test_read_embl_no_header(self):
        "Test reading an embl record without header"
        with open(get_file_path('no_header.embl'), 'rU') as h:
            # plain BioPython reading should fail
            self.assertRaises(ValueError, seqio.read, h)
            h.seek(0)
            # robust reading should work
            record = seqio.read(h, robust=True)
            self.assertEqual("DUMMY.1", record.id)


    def test_read_fasta_no_header(self):
        "Test reading a fasta record without header"
        with open(get_file_path('no_header.fasta'), 'rU') as h:
            # plain BioPython reading should fail
            self.assertRaises(ValueError, seqio.read, h)
            h.seek(0)
            # robust reading should work
            record = seqio.read(h, robust=True)
            self.assertEqual("DUMMY", record.id)


class TestSeqIOZipped(unittest2.TestCase):
    def test_parse_genbank(self):
        "Test parsing a gzipped GenBank file"
        with open(get_file_path('melanin.gbk.gz'), 'rb') as h:
            records = list(seqio.parse(h))
        self.assertEqual(1, len(records))

    def test_read_genbank(self):
        "Test reading a gzipped GenBank file"
        with open(get_file_path('melanin.gbk.gz'), 'rb') as h:
            record = seqio.read(h)
        self.assertEqual("AB070938.1", record.id)

    def test_parse_genbank_path(self):
        "Test parsing a gzipped GenBank file specified by path"
        fname = get_file_path('melanin.gbk.gz')
        records = list(seqio.parse(fname))
        self.assertEqual(1, len(records))

    def test_read_genbank_path(self):
        "Test reading a gzipped GenBank file specified by path"
        fname = get_file_path('melanin.gbk.gz')
        record = seqio.read(fname)
        self.assertEqual("AB070938.1", record.id)
