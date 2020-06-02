from setuptools import setup
from jcast import __version__

setup(
   name='jcast',
   version=__version__,
   description='Jcast retrieves splice junction information and translates into amino acid',
   author='Edward Lau',
   url='https://www.laulab.net',
   author_email='edward.lau@cuanschutz.edu',
   packages=['jcast'],  #same as name
   install_requires=['biopython', 'gtfparse', 'pandas', 'requests', 'tqdm'], #external packages as dependencies
)
