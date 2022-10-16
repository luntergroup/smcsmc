import smcsmc
import filecmp


def test_vcf_to_seg():
    test_vcf = "test/data/test.vcf.gz"
    test_mask = "test/data/test.bed.gz"
    chr = 2

    sample_1 = "ID1"
    sample_2 = "ID2"

    smcsmc.vcf_to_seg(
        [(test_vcf, sample_1), (test_vcf, sample_2)],
        "test/out/test.seg.gz",
        [test_mask, test_mask],
        tmpdir="test/out/tmp",
        key="testconv",
        chroms=[2],
    )
    filecmp.cmp("test/out/test.seg.gz", "test/data/truth/test_conversion.seg.gz")


def test_vcf_to_seg_nomask():
    test_vcf = "test/data/test.vcf.gz"
    test_mask = "test/data/test.bed.gz"
    chr = 2

    sample_1 = "ID1"
    sample_2 = "ID2"

    smcsmc.vcf_to_seg(
        [(test_vcf, sample_1), (test_vcf, sample_2)],
        "test/out/testNoMask.seg.gz",
        tmpdir="test/out/tmp",
        key="testNoMask",
        chroms=[2],
    )
    filecmp.cmp(
        "test/out/testNoMask.seg.gz", "test/data/truth/test_conversion_no_mask.seg.gz"
    )
