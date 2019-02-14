'''
1. Copy a pair input fastq files or a bam file to tests/data/test_input
2. Copy genomo reference fasta as tests/data/reference/NC_000962_3.fasta
3. Copy genomo reference mask array as tests/data/reference/NC_000962_3_repmask.array
4. Copy expected basecall output fasta to tests/data/test_output/expected_output
    python3 tests/test_bwa.py bam  (under CompassCompact, test bam input)
    python3 tests/test_bwa.py fastq  (under CompassCompact, test fastq input)
The test will run the bwa nextflow pipeline and compare the output fasta file with expected fasta file.
'''

import unittest
import subprocess
import os
import shlex
import time
import compare_fasta
import glob
import pathlib
import logging
import sys

class Test_main (unittest.TestCase):
    def setUp(self):
        self.dir = pathlib.Path(os.getcwd())
        self.nf = "main_bwa.nf"
        self.ref = self.dir / "tests/data/reference/NC_000962_3.fasta"
        self.mask_file = self.dir / "tests/data/reference/NC_000962_3_repmask.array"
        self.input_dir = self.dir / "tests/data/test_input/"
        self.pipeline_output = self.dir / "tests/data/test_output/pipeline_output/"
        self.expected_output = self.dir / "tests/data/test_output/expected_output/"

    #@unittest.skip("only skip for debug purpose")
    def test_nf_run_bam(self):
        self.fastq = "false"
        self.pattern = "*.bam"
        self.cmd = ' '.join([
            'nextflow run ', self.nf,
            '--input_dir ', str(self.input_dir) + '/',
            '--output_dir ', str(self.pipeline_output),
            '--pattern ',self.pattern,
            '--ref', str(self.ref),
            '--mask_file', str(self.mask_file),
            '--fastq',self.fastq,
            '-resume', 
            '-profile test_docker'
        ])
        print(self.cmd)
        nextflowprocess = subprocess.Popen(shlex.split(self.cmd))
        nextflowprocess.wait()

    #@unittest.skip("only skip for debug purpose")
    def test_nf_run_fastq(self):
        self.fastq = "true"
        self.pattern = "*_{1,2}.fastq.gz"
        self.cmd = ' '.join([
            'nextflow run ', self.nf,
            '--input_dir ', str(self.input_dir) + '/',
            '--output_dir ', str(self.pipeline_output),
            '--pattern ',self.pattern,
            '--ref', str(self.ref),
            '--mask_file', str(self.mask_file),
            '--fastq',self.fastq,
            '-resume',
            '-profile test_docker'
        ])
        logging.info(self.cmd)
        nextflowprocess = subprocess.Popen(shlex.split(self.cmd))
        nextflowprocess.wait()
        
    def test_nf_run_has_output_folder(self):
        output_dir = self.pipeline_output
        self.assertTrue(os.path.isdir(str(output_dir)))

    def test_nf_run_has_output_files_fasta(self):
        for file in list(self.pipeline_output.glob("*/*consensus.fasta.gz")):
            self.assertTrue(os.path.isfile(str(file)))
    
    def test_nf_run_has_output_files_indel_vcf(self):
        for file in list(self.pipeline_output.glob("*/*basecall_Indel.vcf.gz")):
            self.assertTrue(os.path.isfile(str(file)))

    def test_nf_run_has_output_files_basecall_vcf(self):
        for file in list(self.pipeline_output.glob("*/*basecall.vcf.gz")):
            self.assertTrue(os.path.isfile(str(file)))

    def test_nf_run_has_output_file_basecall_stats(self):
        for file in list(self.pipeline_output.glob("*/*basecallstats.txt")):
            self.assertTrue(os.path.isfile(str(file)))

    def test_nf_run_has_expected_fasta_file(self):
        expected_fasta = []
        pipeline_output_fasta = []
        for file in list(self.pipeline_output.glob("*/*consensus.fasta.gz")):
            expected_fasta.append(str(file))

        for file in list(self.pipeline_output.glob("*/*consensus.fasta.gz")):
            pipeline_output_fasta.append(str(file))
        fasta_as_expected = compare_fasta.fasta_as_expected(expected_fasta[0], pipeline_output_fasta[0])
        self.assertTrue(fasta_as_expected)


if __name__ == '__main__':
    test_type = sys.argv[1]
    if test_type == "bam":
        tests = ['test_nf_run_bam',
                'test_nf_run_has_expected_fasta_file',
                'test_nf_run_has_output_folder',
                'test_nf_run_has_output_files_fasta',
                'test_nf_run_has_output_files_indel_vcf',
                'test_nf_run_has_output_file_basecall_stats',
                'test_nf_run_has_expected_fasta_file']
    if test_type == "fastq":
        tests = ['test_nf_run_fastq',
                'test_nf_run_has_expected_fasta_file',
                'test_nf_run_has_output_folder',
                'test_nf_run_has_output_files_fasta',
                'test_nf_run_has_output_files_indel_vcf',
                'test_nf_run_has_output_file_basecall_stats',
                'test_nf_run_has_expected_fasta_file']
    suite = unittest.TestSuite(map(Test_main, tests))
    unittest.TextTestRunner(verbosity=2).run(suite)
