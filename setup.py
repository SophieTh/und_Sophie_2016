

from setuptools import setup

setup(name='pySRU',
      version='0.5.2',
      description='Python synchrotron undulator calcution',
      author='Sophie Thery',
      author_email='sophie.thery@esrf.fr',
      url='https://github.com/SophieTh/und_Sophie_2016/',
      packages=['pySRU'],
      install_requires=[
                        'numpy',
                        'scipy'
                       ],
      test_suite='tests'
     )



# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# #/*##########################################################################
# #
# # The pySRU library (Synchrotron Radiation Undulator emission in python)
# #
# # Copyright (c) 2016 European Synchrotron Radiation Facility
# #
# # This file is part of pySRU package developed by
# # Sophie Thery, Mark Glass and Manuel Sanchez del Rio
# #
# # Permission is hereby granted, free of charge, to any person obtaining a copy
# # of this software and associated documentation files (the "Software"), to deal
# # in the Software without restriction, including without limitation the rights
# # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# # copies of the Software, and to permit persons to whom the Software is
# # furnished to do so, subject to the following conditions:
# #
# # The above copyright notice and this permission notice shall be included in
# # all copies or substantial portions of the Software.
# #
# # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# # THE SOFTWARE.
# #
# #############################################################################*/
#
# import imp
# import os
# import sys
# import subprocess
#
# # try:
# #     from setuptools import setup
# # except ImportError:
# #     import ez_setup
# #     ez_setup.use_setuptools()
# #     from setuptools import setup
#
# if "bdist_wheel" in sys.argv:
#     from setuptools import setup
#     # from setuptools import Extension
#     # from setuptools.command.build_py import build_py
#     # from distutils.command.install_data import install_data
# else:
#     try:
#         from setuptools import setup
#         # from setuptools import Extension
#         # from setuptools.command.build_py import build_py
#         # from distutils.command.install_data import install_data
#     except ImportError:
#         from distutils.core import setup
#         # from setuptools import Extension
#         # from distutils.sysconfig import get_python_lib
#         # from distutils.command.build_py import build_py
#         # from distutils.command.install_data import install_data
#
#
# NAME = 'pySRU'
#
# VERSION = '1.0.0'
# ISRELEASED = False
#
# DESCRIPTION = 'Synchrotron Radiation Undulator emission'
# README_FILE = os.path.join(os.path.dirname(__file__), 'README.rst')
# LONG_DESCRIPTION = open(README_FILE).read()
# AUTHOR = 'Sophie Thery, Mark Glass and Manuel Sanchez del Rio'
# AUTHOR_EMAIL = 'srio@esrf.eu'
# URL = 'https://github.com/SophieTh/und_Sophie_2016'
# DOWNLOAD_URL = 'https://github.com/SophieTh/und_Sophie_2016'
# MAINTAINER = 'Manuel Sanchez del Rio'
# MAINTAINER_EMAIL = 'srio@esrf.eu'
# LICENSE = 'MIT'
#
# KEYWORDS = (
#     'x-ray'
#     'synchrotron radiation',
#     'emission'
#     'simulation',
# )
#
# CLASSIFIERS = (
#     'Development Status :: 1 - Planning',
#     'Environment :: Console',
#     'Environment :: Plugins',
#     'Programming Language :: Python :: 3',
#     'License :: OSI Approved :: MIT License',
#     'Operating System :: Microsoft :: Windows',
#     'Operating System :: Unix',
#     'Operating System :: MacOS :: MacOS X',
#     'Operating System :: POSIX',
#     'Topic :: Scientific/Engineering :: Physics, Mathematics',
#     'Topic :: Software Development :: Libraries :: Python Modules',
#     'Intended Audience :: Education',
#     'Intended Audience :: Science/Research',
#     'Intended Audience :: Developers',
# )
#
# INSTALL_REQUIRES = (
#     'setuptools',
#     'numpy',
#     'scipy',
# )
#
# SETUP_REQUIRES = (
#     'setuptools',
# )
#
#
# # Return the git revision as a string
# def git_version():
#     """Return the git revision as a string.
#
#     Copied from numpy setup.py
#     """
#     def _minimal_ext_cmd(cmd):
#         # construct minimal environment
#         env = {}
#         for k in ['SYSTEMROOT', 'PATH']:
#             v = os.environ.get(k)
#             if v is not None:
#                 env[k] = v
#         # LANGUAGE is used on win32
#         env['LANGUAGE'] = 'C'
#         env['LANG'] = 'C'
#         env['LC_ALL'] = 'C'
#         out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
#         return out
#
#     try:
#         out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
#         GIT_REVISION = out.strip().decode('ascii')
#     except OSError:
#         GIT_REVISION = "Unknown"
#
#     return GIT_REVISION
#
#
# def write_version_py(filename='pySRU/version.py'):
#     # Copied from numpy setup.py
#     cnt = """
# # THIS FILE IS GENERATED FROM SETUP.PY
# short_version = '%(version)s'
# version = '%(version)s'
# full_version = '%(full_version)s'
# git_revision = '%(git_revision)s'
# release = %(isrelease)s
#
# if not release:
#     version = full_version
#     short_version += ".dev"
# """
#     FULLVERSION = VERSION
#     if os.path.exists('.git'):
#         GIT_REVISION = git_version()
#     elif os.path.exists('pySRU/version.py'):
#         # must be a source distribution, use existing version file
#         version = imp.load_source("pySRU.version", "pySRU/version.py")
#         GIT_REVISION = version.git_revision
#     else:
#         GIT_REVISION = "Unknown"
#
#     if not ISRELEASED:
#         FULLVERSION += '.dev0+' + GIT_REVISION[:7]
#
#     a = open(filename, 'w')
#     try:
#         a.write(cnt % {'version': VERSION,
#                        'full_version': FULLVERSION,
#                        'git_revision': GIT_REVISION,
#                        'isrelease': str(ISRELEASED)})
#     finally:
#         a.close()
#
# #TODO mv tests inside pySRU? pySRU.tests
# PACKAGE_DIR = {""}
#
# PACKAGES = [
#     "pySRU",
#     # "pySRU.tests",
# ]
#
#
# if __name__ == '__main__':
#     write_version_py()
#     setup(
#         name=NAME,
#         version=VERSION,
#         description=DESCRIPTION,
#         long_description=LONG_DESCRIPTION,
#         author=AUTHOR,
#         author_email=AUTHOR_EMAIL,
#         maintainer=MAINTAINER,
#         maintainer_email=MAINTAINER_EMAIL,
#         url=URL,
#         download_url=DOWNLOAD_URL,
#         license=LICENSE,
#         keywords=KEYWORDS,
#         classifiers=CLASSIFIERS,
#         # package_dir=PACKAGE_DIR,
#         packages=PACKAGES,
#         # package_data=PACKAGE_DATA,
#         # extra setuptools args
#         zip_safe=False,  # the package can run out of an .egg file
#         include_package_data=True,
#         install_requires=INSTALL_REQUIRES,
#         setup_requires=SETUP_REQUIRES,
#     )
