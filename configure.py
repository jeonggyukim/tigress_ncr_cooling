#!/usr/bin/env python

import argparse
import re

makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

description = (
    "Prepare custom Makefile and defs.hpp for compiling tigress_ncr_cooling solver"
)
parser = argparse.ArgumentParser(description=description)

cxx_choices = [
    'g++',
    'g++-13',
    'g++-simd',
    'icpx',
    'icpc',
]

def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)
    arg = arg.replace('icx', 'icpx', 1)
    return arg

# --cxx=[name] argument
parser.add_argument(
    '--cxx',
    default='icpx',
    type=c_to_cpp,
    choices=cxx_choices,
    help='select C++ compiler and default set of flags (works w/ or w/o -mpi)')

# -debug argument
parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

# --gcovcmd=[name] argument
parser.add_argument('--gcovcmd',
                    default=None,
                    help='override for command to use to call Gcov utility in Makefile')

# Parse command-line inputs
args = vars(parser.parse_args())

makefile_options = {}
if args['cxx'] == 'g++':
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
elif args['cxx'] == 'g++-13':
    makefile_options['COMPILER_COMMAND'] = 'g++-13'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -I/opt/homebrew/include' # -fopenmp'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
elif args['cxx'] == 'icpx':
    makefile_options['COMPILER_COMMAND'] = 'icpx'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -qopenmp-simd '
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''

# -debug argument
if args['debug']:
    makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g -DATTACH_DEBUGGER'  # -Og
else:
    pass

# --gcovcmd=[name] argument (only modifies Makefile target)
if args['gcovcmd'] is not None:
    makefile_options['GCOV_COMMAND'] = args['gcovcmd']
else:
    makefile_options['GCOV_COMMAND'] = 'gcov'

# Read templates
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()

# Make substitutions
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

# Write output files
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)
