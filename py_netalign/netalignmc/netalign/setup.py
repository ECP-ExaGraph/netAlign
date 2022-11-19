import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

cpp_args = ['-std=c++11', '-fopenmp', '-O3', '-Wall']

ext_modules = [
    Extension(
    'netAlignPY',
        ['netAlignImpl.cpp', 'netAlignPY.cpp'],
        include_dirs=['pybind11/include'],
    language='c++',
    extra_compile_args = cpp_args,
    ),
]

setup(
    name='netAlignPY',
    version='0.0.1',
    author='Arif Khan',
    author_email='arif.khan@pnnl.gov',
    description='Network Alignment',
    ext_modules=ext_modules,
)
