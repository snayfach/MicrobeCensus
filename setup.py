try: from setuptools import setup
except: from distutils.core import setup

setup(
	name = 'MicrobeCensus',
	version = '1.3.0',
	packages = ['microbe_census', 'training', 'tests'],
	package_data={
		'microbe_census': ['data/*', 'bin/*', 'example/*'],
		#'training': ['example/*'],
		'tests': ['data/*']},
	scripts=['scripts/run_microbe_census.py'],
	license = 'GPL',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	url='https://github.com/snayfach/MicrobeCensus',
	install_requires = ['biopython', 'numpy']
)
