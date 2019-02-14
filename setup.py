import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

cpp_args = ['-std=c++14', '-lCGAL', '-lz', '-fopenmp', '-fext-numeric-literals']

ext_modules = [
    Extension(
        'example',
        ['src/example.cpp'],
        include_dirs=['include', 'pybind11/include', "usr/include", "/usr/include/eigen3", "cmake-build-release/_deps/pybind11-src/include"],
        language='c++',
        extra_compile_args = cpp_args,
    ),
]

setup(
    name='example',
    version='0.0.1',
    author='Cliburn Chan',
    author_email='cliburn.chan@duke.edu',
    description='Example',
    ext_modules=ext_modules,
)