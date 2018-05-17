from distutils.core import setup
import numpy 

from Cython.Build import cythonize
from distutils.extension import Extension

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext


class CTypes(Extension): pass

class build_ext(build_ext):
    # https://stackoverflow.com/a/34830639
    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypes)
        return super().build_extension(ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        # remove the python 3 prefixes which we cannot determine
        # at runtime.
        if self._ctypes:
            return ext_name + '.so'
        return super().get_ext_filename(ext_name)

extensions = [
        CTypes("hfof.libhfof", [
            "src/fof.c",  "src/fof64.c", "src/testing.c", "src/periodic.c"],
            include_dirs=["src/", numpy.get_include()],
            extra_compile_args=['-O3', '-g'],
            extra_link_args=['-O3', '-g'],
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

setup(name="hfof", version=find_version("hfof/version.py"),
      author="Peter Creasey",
      author_email="pec27",
      description="",
      url="http://github.com/pec27/hfof",
      zip_safe=False,
      package_dir = {'hfof': 'hfof'},
      packages = [
        'hfof', 'hfof.tests',
      ],
      package_data = {
        'hfof': ['tests/*.dat']
      },
      license='GPL3',
      install_requires=['numpy', 'cython'],
      ext_modules = cythonize(extensions),
      cmdclass={'build_ext' : build_ext}
      )

