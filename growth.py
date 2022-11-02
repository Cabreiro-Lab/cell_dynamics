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
__date__ = 'Nov 2022'
__version__ = '0.3'

from functions.functions import *
import argparse
import pathlib
from scipy.optimize import curve_fit
from scipy.signal import wiener
from time import perf_counter
import time
from functions.functions import *
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
import matplotlib.style as mplstyle
from tqdm import tqdm
import shutil
from multiprocessing import Pool
from multiprocessing import get_context

from scipy import interpolate as ip

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
parser.add_argument('-t',
                    '--threads',
                    default=1,
                    help='Number of threads to use. 1 by default.')
parser.add_argument('-m',
                    '--mode',
                    default='growth',
                    help='Select between the 2 modes: biolog or growth')

args = parser.parse_args()

# this function calculates the AUC of the OD time series
def out_auc_function(file, mode):
    # open the info about df, temperatures, times and OD
    df, temps, OD = biospa_text_opener(os.path.join(ROOT, file))
    OD = OD[3:]

    if check_outliers_temps(temps):
        raise Exception("There are outliers in the temperatures, check your experiment!")
    else:
        pass

    # get time info
    length, time_h, time_span = get_time_h(df)

    ### Fix and interpolate the data
    window = int(round(length / 10, 0))
    # make all the time series start from the same point (0 usually)
    adj_df = df.apply(lambda row: pd.Series(set_series2min(row - np.mean(row[:window]), 0.0), index=time_span), axis=1)

    # smooth the data from adj_df using a wiener filter
    with np.errstate(divide='ignore', invalid='ignore'): # ignore the warnings from wiener filter
        w_filt = adj_df.apply(lambda row: pd.Series(signal.wiener(row, 5), index=time_span), axis=1)
    # fix the NaN rows and substitute them with 0
    w_filt = w_filt.fillna(0)

    growth_rates = w_filt.apply(lambda row: pd.Series(ip.UnivariateSpline(time_h, row, s=0).derivative(1)(time_h), index=time_span),
                    axis=1)

    max_slope = growth_rates.apply(lambda row: pd.Series(row.rolling(window).mean().max(), index=time_span), axis=1).iloc[:,0]

    # calculate AUC from the smoothed data
    auc = w_filt.apply(lambda row: pd.Series(np.trapz(row, time_h), index=time_span), axis=1).iloc[:,0]

    ######## prepare the Output.csv
    # remove RuntimeWarning: divide by zero encountered in log2
    with np.errstate(divide='ignore', invalid='ignore'):
        auc_log2 = np.log2(auc)

    auc_df = pd.DataFrame({f'File': file, f'{OD}_f_AUC': auc, f'{OD}_f_logAUC': auc_log2, f'{OD}_dt_Max': max_slope})
    auc_df.reset_index(inplace=True)
    auc_df = auc_df.rename(columns = {0:'Well'})

    # get the timeseries
    # add metadata to adj_df
    adj_df.insert(0, 'File', file)
    adj_df.insert(1, 'Data', f'{OD}nm')
    adj_df.insert(2, 'Well', adj_df.index)
    # add metadata to w_filt
    w_filt.insert(0, 'File', file)
    w_filt.insert(1, 'Data', f'{OD}nm_f')
    w_filt.insert(2, 'Well', w_filt.index)
    # add metadata to growth_rates
    growth_rates.insert(0, 'File', file)
    growth_rates.insert(1, 'Data', f'{OD}nm_dt')
    growth_rates.insert(2, 'Well', growth_rates.index)

    # concat the dfs
    timeseries_df = pd.concat([adj_df, w_filt, growth_rates], axis=0)

    if mode == 'AUC':
        return auc_df
    elif mode == 'timeseries':
        return timeseries_df


print(f'The folder {os.path.split(args.input)[-1]} will be analysed in the mode \"{args.mode}\".')

# initialise variables
ROOT = args.input
n_threads = int(args.threads)
# general parameters
WAVES = []  # all info about the nm will be stored here to output the most commonly used term for the OD
WAVE = 0  # wave length used in the experiment

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
# create an output folder, overwrite if it exists
if not os.path.exists(os.path.join(args.input, args.output)):
    os.makedirs(os.path.join(args.input, args.output))
else:
    # give option to overwrite or not
    overwrite = input(f'The folder {args.output} already exists. Do you want to overwrite it? (y/n): ')
    if overwrite == 'y':
        shutil.rmtree(os.path.join(args.input, args.output))
        os.makedirs(os.path.join(args.input, args.output))
    else:
        print('Exiting from the run.')
        exit()

# create a folder for the plots
if not os.path.exists(os.path.join(args.input, args.output, 'Plots')):
    os.makedirs(os.path.join(args.input, args.output, 'Plots'))
else:
    # give option to overwrite or not
    overwrite = input(f'The folder {args.output}/Plots already exists. Do you want to overwrite it? (y/n): ')
    if overwrite == 'y':
        shutil.rmtree(os.path.join(args.input, args.output, 'Plots'))
        os.makedirs(os.path.join(args.input, args.output, 'Plots'))
    else:
        print('Plots will be saved in the existing folder.')
        
# LOOP OVER THE FILES
# parallel loop to get the AUCs
with get_context("fork").Pool(10) as p:
    # user starmap to pass multiple arguments to the function
    out_auc_df = pd.concat(list(tqdm(p.starmap(
                                               out_auc_function, 
                                               zip(design.File.to_list(), 
                                               # modify this list to pass different arguments to the function
                                               ['AUC']*len(design.File.to_list()))), 
                                total=len(design.File.to_list()))), axis=0)
p.close()

# parallel loop to get the timeseries
with get_context("fork").Pool(10) as p:
    # user starmap to pass multiple arguments to the function
    out_time_df = pd.concat(list(tqdm(p.starmap(
                                               out_auc_function, 
                                               zip(design.File.to_list(), 
                                               # modify this list to pass different arguments to the function
                                               ['timeseries']*len(design.File.to_list()))), 
                                total=len(design.File.to_list()))), axis=0)
p.close()

# merge the out_auc_df with the design 
out_auc_df = design.merge(out_auc_df, on='File', how='left')
# merge the out_time_df with the design
out_time_df = design.merge(out_time_df, on='File', how='left')

# save the output file in the output folder as a csv file
out_auc_df.to_csv(os.path.join(args.input, args.output, 'Summary.csv'), index=False)
# save the timeseries file in the output folder as a csv file
out_time_df.to_csv(os.path.join(args.input, args.output, 'Timeseries.csv'), index=False)


# to test: python growth.py -i ./test_data/ -t 6

print("All is OK")
