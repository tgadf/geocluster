from distutils.core import setup
import setuptools

setup(
  name = 'geocluster',
  py_modules = ['geocluster', 'plots', 'geoclusterUtils'],
  version = '0.0.1',
  description = 'A Geospatial Clustering Tool',
  long_description = open('README.md').read(),
  author = 'Thomas Gadfort',
  author_email = 'tgadfort@gmail.com',
  license = "MIT",
  url = 'https://github.com/tgadf/geocluster',
  keywords = ['geohash', 'location', 'clustering'],
  classifiers = [
    'Development Status :: 3',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Topic :: Utilities'
  ],
  install_requires=['utils==0.0.1', 'geohash==0.0.1', 'haversine', 'pandas', 'numpy'],
  dependency_links=['git+ssh://git@github.com/tgadf/utils.git#egg=utils-0.0.1', 'git+ssh://git@github.com/tgadf/geohash.git#egg=geohash-0.0.1']
)
 
