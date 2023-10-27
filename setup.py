"""
jastro setup.py file

Inspiration is taken from the donuts package setup
"""
from setuptools import setup

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
    conf = ConfigParser()
    conf.read(['setup.cfg'])
    metadata = dict(conf.items('metadata'))

PACKAGENAME = metadata.get('package_name', 'packagename')
DESCRIPTION = metadata.get('description', '')
AUTHOR = metadata.get('author', '')
AUTHOR_EMAIL = metadata.get('author_email', '')
LICENSE = metadata.get('license', '')
URL = metadata.get('url', 'https://github.com/jmccormac01/jastro')
LONG_DESCRIPTION = open('README.md').read()

from jastro_version import VERSION, RELEASE

# add sphinx build_docs integration
#from sphinx.setup_command import BuildDoc
#cmdclass = {'build_sphinx': BuildDoc}

setup(name=PACKAGENAME,
    version=VERSION,
    description=DESCRIPTION,
    packages=['jastro',],
    install_requires=['numpy>=1.23.5',
                      'pyregion>=2.2.0',
                      'ccdproc>=2.4.0',
                      'sep>=1.2.1',
                      'donuts>=0.3.5',
                      'fitsio>=1.2.0',
                      'astropy>=5.2.1'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov'],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=URL,
    long_description=LONG_DESCRIPTION,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Astronomy',
    ]
)
