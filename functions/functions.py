import os
import pandas as pd
from glob import glob


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


def check_if_equal(list_1, list_2):
    """ Check if both the lists are of same length and if yes then compare
    sorted versions of both the list to check if both of them are equal
    i.e. contain similar elements with same frequency. """
    if len(list_1) != len(list_2):
        return False
    return sorted(list_1) == sorted(list_2)


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

