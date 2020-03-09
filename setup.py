import setuptools

with open('README.md', 'r') as fh:
    long_description=fh.read()

setuptools.setup(name='ramslibs',
      version='0.6.1',
      description='Set of tools for working with RAMS data',
      long_description=long_description,
      long_description_content='text/markdown',
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
