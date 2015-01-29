from setuptools import setup, find_packages

setup(
	name = 'MicrobeCensus',
	version = '1.2.1',
	description = 'A command-line tool for estimating average genome size from shotgun sequence data',
	install_requires=['numpy', 'biopython'],
	license = 'GPL',
	url = 'https://github.com/snayfach/MicrobeCensus',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	classifiers = [
		'Programming Language :: Python',
		'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
		'Development Status :: 5 - Production/Stable',
		'Operating System :: MacOS :: MacOS X',
		'Operating System :: Unix',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],

	keywords='MicrobeCensus metagenomics genomics microbiology ecology',
	scripts=['src/microbe_census.py'],
)
