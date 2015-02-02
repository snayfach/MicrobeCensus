try: from setuptools import setup
except: from distutils.core import setup

setup(
	name = 'test123098',
	version = '1.0.4',
	packages = ['microbe_census'],
	package_data={'microbe_census': ['data/*', 'bin/*']},
	scripts=['scripts/run_microbe_census.py'],
	license = 'GPL',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	url='https://github.com/snayfach/test',
	install_requires = ['biopython', 'numpy']
)
