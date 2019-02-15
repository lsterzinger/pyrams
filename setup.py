from setuptools import setup

setup(name='ramslibs',
      version='0.2',
      description='Set of tools for working with RAMS data',
      author='Lucas Sterzinger',
      author_email='lsterzinger@ucdavis.edu',
      packages=['ramslibs'],
      install_requires=[
          'netCDF4',
          'numpy',
          'matplotlib',
          'metpy',
          'basemap',
          'glob2',
      ]
      )
