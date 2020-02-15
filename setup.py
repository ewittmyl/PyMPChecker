import setuptools
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyMPChecker",
    version="0.1.1",
    author="Chris Klein and Duncan Galloway",
    author_email="duncan.galloway@monash.edu",
    description="A package for querying the Minor Planet database",
    long_description=long_description,
    long_description_content_type="text/markdown",
#    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    package_data={'': ['Monthly_Orbit_Catalogs/*']},
    scripts=glob.glob('scripts/*'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
