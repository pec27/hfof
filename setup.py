from setuptools import setup, Extension


def find_version(path):
    import re
    # path shall be a plain ascii text file.
    s = open(path, 'rt').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              s, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Version not found")

lib = Extension('build.libhfof',
                sources = ["src/fof.c",  "src/fof64.c", "src/testing.c", "src/periodic.c"])

setup(name='hfof', version=find_version("hfof/version.py"),
      author="Peter Creasey",
      author_email="pec27",
      description='Friends of Friends with spatial hashing (hfof)',
      url="http://github.com/pec27/hfof",
      package_dir = {'hfof': 'hfof'},
      packages = ['hfof', 'hfof.tests'],
      license='MIT',
      install_requires=['numpy', 'scipy'],
      ext_modules = [lib],
      test_suite='nose.collector',
      tests_require=['nose'])

