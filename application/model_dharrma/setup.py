### Repris du code de J. CUNHERA dans Solazzi et al., 2021
import subprocess
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize



if not os.path.exists("./lib/"):
   os.makedirs("./lib/")



extensions = [
    Extension("lib.VGfunctions", 
                ["src_Cpp/VGfunctions.pyx"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),

    Extension("lib.RPfunctions", 
                ["src_Cpp/RPfunctions.pyx", "src_Cpp/VGfunctions_src.cpp"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),

    Extension("lib.TTDSPfunctions", 
                ["src_Cpp/TTDSPfunctions.pyx"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),
]

setup(
    ext_modules = cythonize(extensions)
)



subprocess.run(['mv', './src_Cpp/RPfunctions.cpp', './lib/'])
subprocess.run(['mv', './src_Cpp/VGfunctions.cpp', './lib/'])
subprocess.run(['mv', './src_Cpp/TTDSPfunctions.cpp', './lib/'])
