from setuptools import setup


def get_version():
    with open("VERSION", "r") as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as f:
        long_description = f.read()
    return long_description


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]


setup(
    name="metasnek",
    url="https://github.com/beardymcjohnface/metasnek",
    python_requires=">=3.8",
    description="Misc functions for metagenomics pipelines",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Michael Roach",
    author_email="beardymcjohnface@gmail.com",
)
