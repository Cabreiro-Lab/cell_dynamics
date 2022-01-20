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


def biospa_text_opener(file):
    # opens the file
    with open(file, 'r', errors='ignore') as f:
        contents = f.readlines()
    # save the OD
    od = contents[0].split('\n')[0]
    # save the times
    times = contents[2].split('\t')[1:-1]
    # save temperatures
    temps = contents[3].split('\t')[1:-1]
    # save the useful data info
    cosa = contents[4:101]
    # convert it to a pandas object
    df = df_beautify(cosa, times=times)

    return print(f'I could open file {file}')


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
    df.columns = times  # put times as column names

    return df
