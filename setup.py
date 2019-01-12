from os import path
from setuptools import find_packages
import sys
import versioneer
import numpy.distutils.core


f90_files = ['global_module.f90', 'dualpara.f90']
f90_paths = [path.join('pyPIPS', f90_file) for f90_file in f90_files]

ext1 = numpy.distutils.core.Extension(
    name='pyPIPS.dualpara',
    sources=f90_paths,
    extra_f90_compile_args=['-ffree-form', '-ffree-line-length-none']
)


# NOTE: This file must remain Python 2 compatible for the foreseeable future,
# to ensure that we error out properly for people with outdated setuptools
# and/or pip.
min_version = (3, 6)
if sys.version_info < min_version:
    error = """
    pyPIPS does not support Python {0}.{1}.
    Python {2}.{3} and above is required. Check your Python version like so:

    python3 --version

    This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
    Upgrade pip like so:

    pip install --upgrade pip
    """.format(*sys.version_info[:2], *min_version)
    sys.exit(error)
here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]


numpy.distutils.core.setup(
    name='pyPIPS',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=("A collection of python modules and scripts to analyze data from the Purdue/OU"
                 " Portable In situ Precipitation Stations"),
    long_description=readme,
    author="Dan Dawson",
    author_email='meteodan@gmail.com',
    url='https://github.com/Meteodan/pyPIPS',
    packages=find_packages(exclude=['docs', 'tests']),
    entry_points={
        'console_scripts': [
            # 'some.module:some_function',
            ],
        },
    include_package_data=True,
    package_data={
        'pyPIPS': [
            # When adding files here, remember to update MANIFEST.in as well,
            # or else they will not be included in the distribution on PyPI!
            # 'path/to/data_file',
            ]
        },
    install_requires=requirements,
    license="BSD (3-clause)",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    ext_modules=[ext1],
)
