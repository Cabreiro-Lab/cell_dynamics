import os
import pandas as pd
import numpy as np
from glob import glob
from collections import Counter


# functions
def read_my_excel(file_input):
    """
        Reads an Excel from a path with the main design of the experiment.
        It should contain relevant metadata.
        :param file_input: path to an Excel file
        :return: a Pandas object
    """

    xlfile = pd.read_excel(file_input,
                           sheet_name='Design',
                           engine='openpyxl')
    return xlfile


def file_parser(path, pattern='*.txt'):
    """
        Gets the file list within the current folder
        :param path: path where to look for files
        :param pattern: pattern to search
        :return: a list of files
    """
    file_list = []
    for name in map(os.path.basename, glob(os.path.join(path, pattern))):
        file_list.append(name)
    return file_list


def check_if_equal(list_1, list_2):
    """
        Check if both the lists are of same length and if yes then compare
        sorted versions of both the list to check if both of them are equal
        i.e. contain similar elements with same frequency.
    """
    if len(list_1) != len(list_2):
        return False
    return sorted(list_1) == sorted(list_2)


def fix_times(times):
    """
    Takes a vector with the raw times values from the txt and fixes them into a proper shape and rounded values
    :param times: vector of times from the txt file from biospa
    :return: a vector of fixed times
    """
    times = [time_to_sec(val) for val in times if val not in ['Time', '', '0:00:00']]
    length = len(times)
    timestep = round_to(float(times[-1] - times[0]) / (length - 1), 1)
    timemax_min = int((length - 1) * timestep / 60)  # time max in minutes
    time_span = np.linspace(0, timemax_min * 60, length, dtype=np.dtype(int))
    return time_span


def biospa_text_opener(file):
    """
    Main function for opening the txt files from biospa
    :param file: a path to the txt file corresponding to a single plate in your computer
    :return: a tuple of a pandas dataframe, a vector of temperatures, and an OD value for the plate
    """
    # opens the file
    with open(file, 'r', errors='ignore') as f:
        contents = f.readlines()
    # save the OD
    od = contents[0].split('\n')[0]
    # save the times
    times = contents[2].split('\t')[1:-1]
    times = fix_times(times)  # this fixes size and values for the times, rounding them
    # save temperatures
    temps_raw = contents[3].split('\t')[1:-1]
    temps = [float(temp) for temp in temps_raw if temp not in ['0.0', '']]
    temps = np.array(temps)
    # save the useful data info
    temp_df = contents[4:len(contents)]
    # convert it to a pandas object
    df = df_beautify(temp_df, times=times)

    return df, temps, od


def df_beautify(txt_object, times):
    """
    Function to modify pandas dataframes coming from the function biospa_text_opener
    and get them ready to be analysed
    :param times: the times vector from the txt file
    :param txt_object: the txt object that is being read by biospa_txt_opener
    :return: a better pandas dataframe
    """
    df = pd.DataFrame([x.split('\t') for x in txt_object])
    df = df.dropna()
    df = df.set_index(df[0])  # set index as well name
    df = df.drop(0, axis=1)
    df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)  # remove last column as it's a \n
    df = df.replace(r'^\s*$', np.NaN, regex=True)  # replace empty values (from 0:00:00 time values) for NaN
    df.dropna(axis=1, inplace=True)  # remove empty columns
    df.columns = times  # put times as column names
    df = df.apply(pd.to_numeric)  # change type to numeric

    return df


def check_nm(nm):
    """
    Checks the most common wavelength found in the files
    :param nm:
    :return:
    """
    count = Counter(nm)
    if len(count) > 1:
        max_item = list(count.keys())[0]
        item_instances = list(count.values())[0]
        return print(f'Careful! I found more than one wavelength in the files. \n'
                     f'The most common wavelength is {max_item} with {item_instances}')
    else:
        print(f'Wavelength of the experiment is: {nm[0]}')


def round_to(n, precision):
    """
    Function from the original script to round numbers to a desired precision
    (Need to check if I really need this)
    :param n: a float
    :param precision: an integer
    :return:
    """
    # Round a number to desired precision
    correction = 0.5 if n >= 0 else -0.5
    return int(n / precision + correction) * precision


def time_to_sec(time_str):
    """
    Converts a time string from reading the file with biospa_text_opener and converts it to seconds
    :param time_str: a time string like '0:00:30'
    :return: an integer of the total seconds (30 in the case of the example)
    """
    h, m, s = time_str.split(':')
    seconds = int(s) + 60 * int(m) + 3600 * int(h)
    return seconds


def set_series2min(x, thres: float = 0.0):
    """
    This function takes a pandas series and removes all the 0s
    :param x: a Pandas Series object
    :param thres: threshold to use as a minimum value
    :return: returns a numpy array without negative values
    """
    x = x.to_numpy()
    x = np.clip(x, thres, np.inf)
    return x


def growth(x, A, lam, u):
    """
    Parametric logistic growth model.
    Ref: https://www.jstatsoft.org/article/download/v033i07/367
    :param x: series values
    :param A: carrying capacity or max growth
    :param lam: length of lag phase
    :param u: growth rate
    :return: returns the model to be optimised with curve_fit from scipy
    """
    return A / (1 + np.exp((4 * u / A) * (lam - x) + 2))


def check_outliers_temps(temp_vect, thres: float = 0.1):
    """
    Checks if there are outliers in the temperature vector
    :param temp_vect: a numpy array of temperatures
    :param thres: the threshold to detect outliers, in proportion of the average value (0.1 by default)
    :return: a boolean, and the numpy array of the possible outliers
    """
    temps_av = np.average(temp_vect) * thres
    min_temp, max_temp = np.average(temp_vect) - temps_av, np.average(temp_vect) + temps_av
    filter_temps = temp_vect[np.logical_and(temp_vect > max_temp, temp_vect < min_temp)]
    if len(filter_temps) == 0:
        return False
    else:
        return True, filter_temps


def integral(biospa_df, position: str):
    var = biospa_df.loc[position].to_list()

    var = np.array(var)
    var = var - min(var)  # set min values to 0

    return np.trapz(y=var)


# TODO:
#   - introduce functions as dircheck(somedir), to check the existence of folder
#   - start separating functions by classes (reading files, converting, analysing...)
#   that will help to create different files and make the code go faster depending on the options
