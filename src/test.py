import unittest
import tempfile
import ags_functions
import shutil
import os

class ReadList(unittest.TestCase):
	""" check whether file is read into list properly """
	
	def setUp(self): # create test files in tmp directory
		self.exp_values = range(10)
		self.dir = tempfile.mkdtemp()
		self.inpath = os.path.join(self.dir, 'tmp.txt')
		for i in range(10):
			open(self.inpath, 'a').write(str(i)+'\n')
		
	def test_read_list(self):
		obs_values = ags_functions.read_list(self.inpath, header=False, dtype='int')
		for exp, obs in zip(self.exp_values, obs_values):
			self.assertEqual(exp, obs)

	def tearDown(self): # clean up tmp directory
		shutil.rmtree(self.dir)


class FileType(unittest.TestCase):
	""" check whether filetype is correctly determined """

	def setUp(self): # create test files in tmp directory
		self.dir = tempfile.mkdtemp()
		self.values = [['file.fq', 'fastq',
		               ('@HWUSI-EAS574_102539073:1:100:10000:12882/1',
						'AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC',
						'+HWUSI-EAS574_102539073:1:100:10000:12882/1',
						'GGGGFGFFGGAGDFG=EDEEDEBEEEEEEEEEEEAB?B?BEEE')],
					   ['file.fa', 'fasta',
		               ('>HWUSI-EAS574_102539073:1:100:10000:12882/1',
						'AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC')],
					   ['file.txt', None,
		               ('some random text file')]]
		for name, type, values in self.values:
			inpath = os.path.join(self.dir, name)
			infile = open(inpath, 'w')
			for value in values: infile.write(value+'\n')

	def test_detect_filetype(self):
		# fastq
		inpath = os.path.join(self.dir, self.values[0][0])
		type = ags_functions.auto_detect_file_type(inpath)
		self.assertEqual(type, self.values[0][1])
		# fasta
		inpath = os.path.join(self.dir, self.values[1][0])
		type = ags_functions.auto_detect_file_type(inpath)
		self.assertEqual(type, self.values[1][1])
		# neither
		inpath = os.path.join(self.dir, self.values[2][0])
		with self.assertRaises(SystemExit):
			ags_functions.auto_detect_file_type(inpath)

	def tearDown(self): # clean up tmp directory
		shutil.rmtree(self.dir)


class QualEncode(unittest.TestCase):
	""" check whether fastq quality encoding is correctly determined """

	def setUp(self): # create test files in tmp directory
		self.dir = tempfile.mkdtemp()
		self.values = [['file.sanger', 'sanger',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIII""")],
					   ['file.solexa', 'solexa',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATCCCC""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						""";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""")],
					   ['file.illumina', 'illumina',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCA""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""")]]
		for name, type, values in self.values:
			inpath = os.path.join(self.dir, name)
			infile = open(inpath, 'w')
			for value in values: infile.write(value+'\n')

	def test_detect_qualtype(self):
		for filename, exptype, seqs in self.values:
			inpath = os.path.join(self.dir, filename)
			obstype = ags_functions.auto_detect_fastq_format(inpath)
			self.assertEqual(obstype, exptype)

	def tearDown(self): # clean up tmp directory
		shutil.rmtree(self.dir)


class QualEncode(unittest.TestCase):
	""" check whether fastq quality encoding is correctly determined """

	def setUp(self): # create test files in tmp directory
		self.dir = tempfile.mkdtemp()
		self.values = [['file.sanger', 'sanger',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIII""")],
					   ['file.solexa', 'solexa',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATCCCC""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						""";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""")],
					   ['file.illumina', 'illumina',
		               ("""@HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCA""",
						"""+HWUSI-EAS574_102539073:1:100:10000:12882/1""",
						"""@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""")]]
		for name, type, values in self.values:
			inpath = os.path.join(self.dir, name)
			infile = open(inpath, 'w')
			for value in values: infile.write(value+'\n')

	def test_detect_qualtype(self):
		for filename, exptype, seqs in self.values:
			inpath = os.path.join(self.dir, filename)
			obstype = ags_functions.auto_detect_fastq_format(inpath)
			self.assertEqual(obstype, exptype)

	def tearDown(self): # clean up tmp directory
		shutil.rmtree(self.dir)





if __name__ == '__main__':
	unittest.main()