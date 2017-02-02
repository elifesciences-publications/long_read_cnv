from setuptools import setup, find_packages

import os 

def locate_packages():
    packages = ['long_read_pipeline']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages

setup(
    name="lrcnv",
    version="1.0",
    packages=locate_packages(),
    author="James Boocock",
    author_email="james.boocock@ucla.edu",
    description="Causal Effects from Summary Statstics (CESS)",
    license="Mit",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'lr_cnv = long_read_pipeline.main:main'
        ]
        },
    url="github.com/theboocock/lrcnv",
    use_2to3=True
)

