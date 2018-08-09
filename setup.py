from setuptools import setup, Extension

import numpy
from distutils.command.build_ext import build_ext as build_ext_base

class CTypes(Extension): pass

class build_ext(build_ext_base):
    # based on https://stackoverflow.com/a/34830639
    # modified base class invocation to support python 2
    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypes)
        return build_ext_base.build_extension(self, ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return build_ext_base.get_export_symbols(self, ext)

    def get_ext_filename(self, ext_name):
        # remove the python 3 prefixes which we cannot determine
        # at runtime, and also remove python2 ext name prefixes
        if self._ctypes:
            return ext_name.split('.')[-1] + '.so'
        return build_ext_base.get_ext_filename(self, ext_name)

extensions = [
        CTypes("hfof.libhfof", [
            "src/fof.c",  "src/fof64.c", "src/testing.c", "src/periodic.c"],
            include_dirs=["src/", numpy.get_include()],
            extra_compile_args=['-O2', '-g', '-std=c99'],
            extra_link_args=['-O2', '-g', '-std=c99'],
            libraries = ['m'],
            )
        ]

def find_version(path):
    import re
    # path shall be a plain ascii text file.
    s = open(path, 'rt').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              s, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Version not found")

setup(name='hfof', version=find_version("hfof/version.py"),
      author="Peter Creasey",
      author_email="pec27",
      description='Friends of Friends with spatial hashing (hfof)',
      url="http://github.com/pec27/hfof",
      zip_safe=False,
      package_dir = {'hfof': 'hfof'},
      packages = ['hfof', 'hfof.tests'],
      package_data = {
        'hfof': ['tests/*.dat']
      },
      license='MIT',
      install_requires=['numpy', 'scipy'],
      ext_modules = extensions,
      test_suite='nose.collector',
      tests_require=['nose'])

