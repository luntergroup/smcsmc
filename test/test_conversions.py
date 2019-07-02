import smcsmc
import pdb

def test_vcf_to_seg():
    test_vcf = "data/chr{}.vcf.gz"
    test_mask = "data/mask.chr{}.txt.gz"
    chr = 2

    sample_1 = "S_Yoruba-1"
    sample_2 = "S_Mbuti-1"

    smcsmc.vcf_to_seg(
            [(test_vcf.format(chr), sample_1),
             (test_vcf.format(chr), sample_2)],
            [test_mask, 
             test_mask],
            "out/test.seg.gz",
            tmpdir = "out/tmp",
            key = "testconv",
            chroms = [2]) 

test_vcf_to_seg()
