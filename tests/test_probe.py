from src.probe import *


class TestInterval:
    def test_str_nullIntervalReturnsEmptyString(self):
        interval = ProbeInterval()

        actual = str(interval)
        expected = ""

        assert actual == expected

    def test_str_nonNullIntervalReturnsExpectedString(self):
        interval = ProbeInterval(1, 5)

        actual = str(interval)
        expected = "[1,5)"

        assert actual == expected

    def test_len_isZeroLengthReturnsZero(self):
        interval = ProbeInterval(2, 2)

        actual = len(interval)
        expected = 0

        assert actual == expected

    def test_len_isTwoReturnsTwo(self):
        interval = ProbeInterval(4, 6)

        actual = len(interval)
        expected = 2

        assert actual == expected

    def test_len_isNullReturnsZero(self):
        interval = ProbeInterval()

        actual = len(interval)
        expected = 0

        assert actual == expected

    def test_isNull_nullIntervalReturnsTrue(self):
        assert ProbeInterval(-1, -1).is_null()

    def test_isNull_nonNullIntervalReturnsFalse(self):
        assert not ProbeInterval(0, -1).is_null()

    def test_fromString_emptyStringReturnsNone(self):
        string = ""

        actual = ProbeInterval.from_string(string)

        assert actual is None

    def test_fromString_invalidStringReturnsNone(self):
        string = "[10,20]"

        actual = ProbeInterval.from_string(string)

        assert actual is None

    def test_fromString_validStringWithSpaceAfterCommaReturnsNone(self):
        string = "[10, 20)"

        actual = ProbeInterval.from_string(string)

        assert actual is None

    def test_fromString_validStringReturnsInterval(self):
        string = "[10,20)"

        actual = ProbeInterval.from_string(string)
        expected = ProbeInterval(10, 20)

        assert actual == expected


class TestProbeHeader:
    def test_equality_equalReturnsTrue(self):
        p1 = ProbeHeader(sample="foo")
        p2 = ProbeHeader(sample="foo")

        assert p1 == p2

    def test_equality_notEqualReturnsFalse(self):
        p1 = ProbeHeader(interval=ProbeInterval(2, 5))
        p2 = ProbeHeader(interval=ProbeInterval(2, 4))

        assert p1 != p2

    def test_fromString_emptyStringReturnsEmptyProbeHeader(self):
        string = ""

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader()

        assert actual == expected

    def test_fromString_allFieldsInStringReturnsProbeHeaderWithAllFields(self):
        string = ">CHROM=1;SAMPLE=CFT073;POS=1;REF_LEN=25;IV=[0,72);SVTYPE=INDEL;MEAN_FWD_COVG=2;MEAN_REV_COVG=3;GT_CONF=10.9922;" \
                 "COVERAGE=13;PVID=42;NB_ALL=3;ALL_ID=1;" \
                 "NB_DIFF_ALL_SEQ=10;ALL_SEQ_ID=4;NB_OF_SAMPLES=7;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            pos=1,
            ref_length=25,
            interval=ProbeInterval(0, 72),
            svtype="INDEL",
            gt_conf=10.9922,
            coverage=13,
            pangenome_variation_id=42,
            number_of_alleles=3,
            allele_id=1,
            number_of_different_allele_sequences=10,
            allele_sequence_id=4,
            nb_of_samples=7,
        )

        assert actual == expected

    def test_fromString_someFieldsInStringReturnsProbeHeaderWithSomeFields(self):
        string = ">CHROM=1;SAMPLE=CFT073;SVTYPE=INDEL;MEAN_FWD_COVG=2;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            svtype="INDEL",
            gt_conf=10.9922,
        )

        assert actual == expected

    def test_fromString_stringWithInvalidFieldReturnsProbeHeaderWithOnlyValidFields(
            self
    ):
        string = ">CHROM=1;SAMPLE=CFT073;SVTYPE=INDEL;MEAN_FWD_COVG=2;INVALID=foo;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            svtype="INDEL",
            gt_conf=10.9922,
        )

        assert actual == expected

    def test_str_emptyProbeHeaderReturnsEmptyString(self):
        header = ProbeHeader()

        actual = str(header)
        expected = ""

        assert actual == expected

    def test_str_singleVariableProbeHeaderReturnsStringWithOneField(self):
        header = ProbeHeader(svtype="SNP")

        actual = str(header)
        expected = ">SVTYPE=SNP;"

        assert actual == expected

    def test_str_allFieldsInProbeHeaderReturnsStringWithAllFields(self):
        header = ProbeHeader(
            sample="CFT073",
            chrom="1",
            pos=1,
            ref_length=25,
            interval=ProbeInterval(0, 72),
            svtype="INDEL",
            gt_conf=10.9922,
            coverage=13,
            pangenome_variation_id=42,
            number_of_alleles=3,
            allele_id=1,
            number_of_different_allele_sequences=10,
            allele_sequence_id=4,
            nb_of_samples=7,
        )

        actual = str(header)
        expected = ">CHROM=1;SAMPLE=CFT073;POS=1;REF_LEN=25;IV=[0,72);SVTYPE=INDEL;GT_CONF=10.9922;" \
                   "COVERAGE=13;PVID=42;NB_ALL=3;ALL_ID=1;" \
                   "NB_DIFF_ALL_SEQ=10;ALL_SEQ_ID=4;NB_OF_SAMPLES=7;"

        assert actual == expected

        assert ProbeHeader.from_string(actual) == header


class TestProbe:
    def test_str_emptyProbeReturnsEmptyString(self):
        probe = Probe()

        actual = str(probe)
        expected = ""

        assert actual == expected

    def test_str_emptyFullSequenceReturnsHeaderWithNewline(self):
        probe = Probe(ProbeHeader(chrom="3"))

        actual = str(probe)
        expected = ">CHROM=3;\n"

        assert actual == expected

    def test_str_fullProbeReturnsHeaderAndSequence(self):
        probe = Probe(ProbeHeader(pos=4, chrom="3"), full_sequence="foo")

        actual = str(probe)
        expected = ">CHROM=3;POS=4;\nfoo"

        assert actual == expected

        assert Probe.from_string(actual) == probe

    def test_equality_twoEqualProbesReturnsTrue(self):
        p1 = Probe(
            ProbeHeader(sample="foo", interval=ProbeInterval(1, 2)), full_sequence="bar"
        )
        p2 = Probe(
            ProbeHeader(sample="foo", interval=ProbeInterval(1, 2)), full_sequence="bar"
        )

        assert p1 == p2

    def test_equality_twoNonEqualProbesReturnsFalse(self):
        p1 = Probe(
            ProbeHeader(sample="foo", interval=ProbeInterval(1, 2)),
            full_sequence="barr",
        )
        p2 = Probe(
            ProbeHeader(sample="foo", interval=ProbeInterval(1, 2)), full_sequence="bar"
        )

        assert p1 != p2

    def test_getLeftFlank_emptyFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=ProbeInterval(4, 5))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = ""

        assert actual == expected

    def test_getLeftFlank_noLeftFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=ProbeInterval(0, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = ""

        assert actual == expected

    def test_getLeftFlank_singleBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=ProbeInterval(1, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = "a"

        assert actual == expected

    def test_getLeftFlank_multiBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=ProbeInterval(3, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = "abc"

        assert actual == expected

    def test_getRightFlank_emptyFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=ProbeInterval(4, 5))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = ""

        assert actual == expected

    def test_getRightFlank_noRightFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=ProbeInterval(2, 7))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = ""

        assert actual == expected

    def test_getRightFlank_singleBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=ProbeInterval(1, 6))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = "g"

        assert actual == expected

    def test_getRightFlank_multiBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=ProbeInterval(3, 4))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = "efg"

        assert actual == expected

    def test_getCoreSequence_emptyFullSequenceReturnsEmptyString(self):
        header = ProbeHeader(interval=ProbeInterval(3, 4))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_getCoreSequence_intervalForDeletionReturnsEmptyString(self):
        header = ProbeHeader(interval=ProbeInterval(3, 3))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseReturnsSingleBase(self):
        header = ProbeHeader(interval=ProbeInterval(3, 4))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "d"

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseAtStartOfSequenceReturnsSingleBase(
            self
    ):
        header = ProbeHeader(interval=ProbeInterval(0, 1))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "a"

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseAtEndOfSequenceReturnsSingleBase(
            self
    ):
        header = ProbeHeader(interval=ProbeInterval(6, 7))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "g"

        assert actual == expected

    def test_getCoreSequence_intervalForMultiBaseSequenceReturnsMultiBase(self):
        header = ProbeHeader(interval=ProbeInterval(2, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "cde"

        assert actual == expected

    def test_getCoreSequence_intervalOutOfRangeReturnsEmptyString(self):
        header = ProbeHeader(interval=ProbeInterval(8, 9))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_fromString_emptyStringReturnsEmptyProbe(self):
        string = ""

        actual = Probe.from_string(string)
        expected = Probe()

        assert actual == expected

    def test_fromString_headerOnlyStringReturnsProbeWithNoFullSequence(self):
        string = ">CHROM=1;IV=[3,5);"

        actual = Probe.from_string(string)
        expected = Probe(header=ProbeHeader(chrom="1", interval=ProbeInterval(3, 5)))

        assert actual == expected

    def test_fromString_headerAndEmptySequenceInStringReturnsProbeWithNoFullSequence(
            self
    ):
        string = ">CHROM=1;IV=[3,5);\n"

        actual = Probe.from_string(string)
        expected = Probe(header=ProbeHeader(chrom="1", interval=ProbeInterval(3, 5)))

        assert actual == expected

    def test_fromString_headerAndSequenceInStringReturnsFullProbe(self):
        string = ">CHROM=1;IV=[3,5);\nfoo"

        actual = Probe.from_string(string)
        expected = Probe(
            header=ProbeHeader(chrom="1", interval=ProbeInterval(3, 5)),
            full_sequence="foo",
        )

        assert actual == expected

    def test_gtConf(self):
        probe = Probe(header=ProbeHeader(gt_conf=5.5))

        actual = probe.gt_conf
        expected = 5.5

        assert actual == expected

    def test_chrom(self):
        probe = Probe(header=ProbeHeader(chrom="chrom1"))

        actual = probe.chrom
        expected = "chrom1"

        assert actual == expected

    def test_pos(self):
        probe = Probe(header=ProbeHeader(pos=7))

        actual = probe.pos
        expected = 7

        assert actual == expected

    def test_interval(self):
        probe = Probe(header=ProbeHeader(interval=ProbeInterval(4, 6)))

        actual = probe.interval
        expected = ProbeInterval(4, 6)

        assert actual == expected

    def test___get_interval_or_default_interval_if_none___valid_interval(self):
        probe = Probe(header=ProbeHeader(interval=ProbeInterval(4, 6)))

        actual = probe.get_interval_or_default_interval_if_none()
        expected = ProbeInterval(4, 6)

        assert actual == expected

    def test___get_interval_or_default_interval_if_none___none_interval(self):
        probe = Probe()

        actual = probe.get_interval_or_default_interval_if_none()
        expected = ProbeInterval()

        assert actual == expected
