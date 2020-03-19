import setuptools
import re
VERSIONFILE = "ramslibs/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

with open('README.md', 'r') as fh:
    long_description=fh.read()

setuptools.setup(name='ramslibs',
      version=verstr,
      description='Set of tools for working with RAMS data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/lsterzinger/ramslibs',
      author='Lucas Sterzinger',
      author_email='lsterzinger@ucdavis.edu',
      packages=setuptools.find_packages(),
      classifiers= [
          'License :: OSI Approved :: MIT License',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
          'Development Status :: 4 - Beta'
      ])
