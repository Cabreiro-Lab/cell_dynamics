#
#     Copyright (C) 2022, Daniel Martínez Martínez
#

"""
Python script to transform the output data from the Biospa machine into
a format that can be treated with any other programming language.
The main idea comes from the script Growth.py from Povilas repository,
but this implementation is entirely mine.
"""

__author__ = 'Daniel Martínez Martínez'
__copyright__ = 'Copyright (C) 2022 Daniel Martínez Martínez'
__license__ = 'GNU Affero General Public License Version 3'
__email__ = 'dmartimarti **AT** gmail.com'
__maintainer__ = 'Daniel Martínez Martínez'
__status__ = 'Test'
__date__ = 'Jan 2022'
__version__ = '0.1'

from functions.functions import *
import argparse
import pathlib

# start parser
parser = argparse.ArgumentParser()

# define arguments to pass to the script
# here comes one difference, I'll read a path instead of a Design file
parser.add_argument('-i',
                    '--input',
                    type=lambda p: pathlib.Path(p).absolute(),
                    help='Input path where the design file is located (in .xlsx format)',
                    required=True)
parser.add_argument('-f',
                    '--file',
                    default='Design.xlsx',
                    help='Name of the input file. "Design.xlsx" by default.')
parser.add_argument('-o',
                    '--output',
                    default='Output',
                    help='Output folder to save the analysis')
parser.add_argument('-m',
                    '--mode',
                    default='growth',
                    help='Select between the 2 modes: biolog or growth')

args = parser.parse_args()

print(f'The folder {os.path.split(args.input)[-1]} will be analysed in the mode \"{args.mode}\".')

# initialise variables
ROOT = args.input

# read Excel file with Pandas

design = read_my_excel(os.path.join(args.input, args.file))

# first, check if the files from the file and your computer are the same

des_files = design['File'].to_list()
files_in_system = file_parser(path=args.input,
                              pattern='*.txt')

if check_if_equal(des_files, files_in_system):
    pass
else:
    raise Exception("The files in the folder are not the same as "
                    "the files within your design file!\n"
                    "Exiting from the run.")


# from this point, I need to open the files, calculate their stuff, plot them, and save the relevant
# info in an output csv file
for file in files_in_system:
    biospa_text_opener(os.path.join(ROOT, file))

print("All is OK")
