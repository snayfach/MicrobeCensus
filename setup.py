from distutils.core import setup
from distutils.command import install_lib
from distutils import log
import os

class my_install_lib(install_lib.install_lib):
	""" Make all packaged binaries executable """
	def run(self):
		install_lib.install_lib.run(self)
		for fn in self.get_outputs():
			if fn.split('/')[-2] == 'bin':
				mode = ((os.stat(fn).st_mode) | 0555) & 07777
				log.info("changing mode of %s to %o", fn, mode)
				os.chmod(fn, mode)

setup(
	name = 'MicrobeCensus',
	version = '1.0.6',
	description = 'Estimation of average genome size from metagenomic data',
	packages = ['microbe_census', 'tests'],
	package_data={
		'microbe_census': ['data/*', 'bin/*', 'example/*'],
		'tests': ['data/*']},
	scripts=['scripts/run_microbe_census.py'],
	license = 'GPL',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	url='https://github.com/snayfach/MicrobeCensus',
	install_requires = ['biopython', 'numpy'],
	cmdclass={'install_lib':my_install_lib}
)

