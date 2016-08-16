from setuptools import setup

setup(name='pySRU',
      version='0.5',
      description='Python synchrotron undulator calcution',
      author='Sophie Thery',
      author_email='sophie.thery@esrf.fr',
      url='https://github.com/SophieTh/und_Sophie_2016/',
      packages=['pySRU'],
      install_requires=[
                        'numpy',
                        'scipy'
                       ]
     )
