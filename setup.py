""" setup module """

from setuptools import setup, find_packages
import os.path

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='jcast',
    version="0.3.2",
    description='JCAST retrieves splice junction information from RNA-seq dat aand translates amino acids sequences',

    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/ed-lau/jcast',

    author='Edward Lau, Maggie Lam',
    author_email='edward.lau@cuanschutz.edu',

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='scientific multi-omics isoform mass-spectrometry',  # Optional

    packages=find_packages(),

    python_requires='>=3.6, <4',

    install_requires=['biopython>=1.78,<2',
                      'gtfparse>=1,<2',
                      'pandas>=1.0,<2',
                      'requests>=2,<3',
                      'tqdm>=4,<5',
                      'scikit-learn>=0.24,<0.3',
                      'matplotlib>=3.4,<4',
                      ],  # external packages as dependencies
    entry_points={
        'console_scripts': ['jcast=jcast.main:main',
                            ],
    },

    project_urls={
        'Source': 'https://github.com/ed-lau/jcast',
        'Maggie Lam Lab': 'http://www.maggielab.org',
    },

    data_files=[('tests',
                 [os.path.join('tests', 'data', 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'),
                  os.path.join('tests', 'data', 'genome', 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'),
                  os.path.join('tests', 'data', 'rmats', 'A3SS.MATS.JC.txt'),
                  os.path.join('tests', 'data', 'rmats', 'A5SS.MATS.JC.txt'),
                  os.path.join('tests', 'data', 'rmats', 'RI.MATS.JC.txt'),
                  os.path.join('tests', 'data', 'rmats', 'SE.MATS.JC.txt'),
                  os.path.join('tests', 'data', 'rmats', 'MXE.MATS.JC.txt'),
                  ]),
                ],
)
