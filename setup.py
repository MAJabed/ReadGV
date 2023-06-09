from setuptools import find_packages, setup
import io 
# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='readgv',
    packages=find_packages(),
    version='0.1.1',
    author="Mohammed A. Jabed,",
    author_email="jabed.abu@gmail.com",
    keywords=['vasp', 'gaussian', 'compchem', 'dft','tddft','electronic structure'],
    install_requires=['numpy', 'regex'], 
    python_requires='>=3.0',
    description='A python program to parse and post-processing of Gaussian and VASP outputs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'readgv = readgv.readgv:main',
        ],
    }
)
