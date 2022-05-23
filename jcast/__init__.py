# -*- coding: utf-8 -*-

""" init """
from os.path import dirname, basename, isfile, join
import glob
modules = glob.glob(join(dirname(__file__), "*.py"))
modules = [ basename(f)[:-3] for f in modules if isfile(f) 
			and not f.endswith("__init__.py") 
			and not f.endswith("setup.py") 
			and not f.endswith("__main__.py")]
__all__ = modules

__version_info__ = ('0', '3', '4')
__version__ = '.'.join(__version_info__)
