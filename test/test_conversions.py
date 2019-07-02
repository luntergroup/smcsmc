import smcsmc
import pdb
import filecmp

def test_vcf_to_seg():
    test_vcf = "data/test.vcf.gz"
    test_mask = "data/test.bed.gz"
    chr = 2

    sample_1 = "ID1"
    sample_2 = "ID2"

    smcsmc.vcf_to_seg(
            [(test_vcf, sample_1),
             (test_vcf, sample_2)],
            "out/test.seg.gz",
            [test_mask, 
             test_mask],
            tmpdir = "out/tmp",
            key = "testconv",
            chroms = [2]) 
    filecmp.cmp('out/test.seg.gz', 'data/truth/test_conversion.seg.gz')

def test_vcf_to_seg_nomask():
    test_vcf = "data/test.vcf.gz"
    test_mask = "data/test.bed.gz"
    chr = 2

    sample_1 = "ID1"
    sample_2 = "ID2"

    smcsmc.vcf_to_seg(
            [(test_vcf, sample_1),
             (test_vcf, sample_2)],
            "out/testNoMask.seg.gz",
            tmpdir = "out/tmp",
            key = "testNoMask",
            chroms = [2]) 
    filecmp.cmp('out/testNoMask.seg.gz', 'data/truth/test_conversion_no_mask.seg.gz')


test_vcf_to_seg()
#test_vcf_to_seg_nomask()
