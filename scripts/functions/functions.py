import os
import pandas as pd
import numpy as np
from glob import glob
from collections import Counter
import matplotlib.style as mplstyle
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from itertools import product
import seaborn as sns


# functions
def read_design(file_input):
    """
        Reads an Excel from a path with the main design of the experiment.
        It should contain relevant metadata.

        Parameters 
        ----------
        file_input : str
            Path to the Excel file with the main design of the experiment.

        Returns
        -------
        df : pandas.DataFrame
            A pandas DataFrame with the main design of the experiment.
    """

    xlfile = pd.read_excel(file_input,
                           sheet_name='Design',
                           engine='openpyxl')
    return xlfile


def get_sheet_names(file):
    """
    Function that gets the sheet names from an Excel file

    Parameters
    ----------
    file : str
        Path to the Excel file with the main design of the experiment.

    Returns
    -------
    sheet_names : list
        A list with the names of the sheets in the Excel file.
    """
    xlfile = pd.ExcelFile(file)
    # sheet names if they don't start by '_'
    sheet_names = [sheet for sheet in xlfile.sheet_names if sheet[0] != '_']
    return sheet_names

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
        Check if the elements contained in the list_1 are present in the list_2. If true, returns True. If false, returns False.
    """
    return set(list_1).issubset(set(list_2))


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

def get_well_names(num_letters, num_numbers):
    """
    Function that returns a list of all possible combinations of letters (from A to a specified number) and numbers (from 1 to a specified number)
    Parameters
    ----------
    num_letters : int
        The number of letters to be used
    num_numbers : int
        The number of numbers to be used
    Returns
    -------
    list
        A list of all possible combinations of letters (from A to a specified number) and numbers (from 1 to a specified number)
    """
    # get the letters
    letters = [chr(i) for i in range(65, 65+num_letters)]
    # get the numbers
    numbers = [str(i) for i in range(1, num_numbers+1)]
    # get all the combinations of letters and numbers
    well_names = [i+j for i, j in product(letters, numbers)]
    return well_names

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
    
    # get the index of the line where the times start
    for i, line in enumerate(contents):
        if line[:4] == 'Time':
            time_line = i
            break
    
    # save the times
    times = contents[time_line].split('\t')[1:-1]
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
    df = df.set_index(df[0])  # set index as well name
    df = df.drop(0, axis=1)
    df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)  # remove last column as it's a \n
    df = df.dropna()  # remove empty rows
    df = df.replace(r'^\s*$', np.NaN, regex=True)  # replace empty values (from 0:00:00 time values) for NaN
    df.dropna(axis=1, inplace=True)  # remove empty columns
    df.columns = times  # put times as column names
    df = df.apply(pd.to_numeric)  # change type to numeric

    # check if we have missing wells in the dataframe and add them if needed
    # set of letters from A to H
    # get the index and separate letters and numbers, save the unique instances
    df_index = df.index
    letters = [i[0] for i in df_index]
    numbers = [i[1:] for i in df_index]
    # sort the list of letters and numbers
    letters = sorted(list(set(letters)))
    numbers = sorted(list(set(numbers)))
    combinations = get_well_names(len(letters), len(numbers))

    # check if there is a missing well, and if so, add it with 0 values
    if len(df.index) != len(combinations):
        missing_wells = list(set(combinations) - set(df.index))
        for well in missing_wells:
            df.loc[well] = 0
    else:
        pass

    # sort the dataframe by index following the order of the combinations
    df = df.reindex(combinations)

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
    This function takes a pandas series and clip the values between 0 and Inf
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


def get_time_h(df):
    # time_t is the time from the experiment in format 0:00:00
    time_t = df.columns.to_list()
    time_t = [t+52 for t in time_t]
    # length is the number of time points
    length = len(time_t)
    # timestep is the time between each time point
    timestep = round_to(float(time_t[-1] - time_t[0]) / (length - 1), 1)
    # timemax_min is the total time in minutes
    timemax_min = int((length - 1) * timestep / 60)  # time max in mins
    # timemax_remmin is the remainder of minutes
    timemax_h, timemax_remmin = divmod(timemax_min, 60)
    # time_span is the time in seconds
    time_span = np.linspace(0, timemax_min * 60, length, dtype=np.dtype(int))
    # time_h is the time in hours
    time_h = time_span / 3600.0
    return length, time_h, time_span


# def plot_individual_plate(data, title, out_name, time_h, save=False):
#     """
#     Plots the data as a grid of 8x12. It takes a Pandas dataframe as input with Wells as index and time as columns.

#     :param data: data to plot
#     :param title: title of the plot
#     :param out_name: name of the output file
#     :param save: if True, saves the plot as a pdf file
#     """

#     # get index from data and separate it into numbers and letters, save only a unique list of both
#     plt.ioff()
#     mplstyle.use('fast')
#     index = data.index.to_list()
#     letters = list(set([i[0] for i in index]))
#     numbers =list(set([int(i[1:]) for i in index]))

#     # sort the lists
#     letters.sort()
#     numbers.sort()

#     let_len = len(letters)
#     num_len = len(numbers)

#     # plot data as an 8x12 grid, all plots share the same y and x axis, place a unique axis for each row and column
#     # print numbers as the general column name and letters as the general row name
#     fig, axs = plt.subplots(let_len, num_len, sharex=True, sharey=True, figsize=(num_len, let_len))
#     fig.suptitle(title)
#     for ax, col in zip(axs[0], numbers):
#         ax.set_title(col)
#     for ax, row in zip(axs[:, len(numbers)-1], letters):
#         ax.set_ylabel(row, rotation=0, size='large', labelpad=10)
#         ax.yaxis.set_label_position("right")
#     for i, ax in enumerate(axs.flat):
#         ax.plot(time_h, data.iloc[i])
#         ax.set_xticks(np.arange(0, 24, 12))
#         ax.set_yticks(np.arange(0, max(data.max(axis=1)), 0.3))
    
#     if save:
#         # save the plot in pdf format
#         plt.savefig(f'{out_name}.pdf', format='pdf')


def plot_individual_plate_plotly(data, title, out_name, time_h, save=False):
    """
    Plots the data as a grid of 8x12 using plotly. It takes a Pandas dataframe as input with Wells as index and time as columns.

    Parameters
    ----------
    data : Pandas dataframe
        data to plot
    title : str
        title of the plot
    out_name : str
        name of the output file
    time_h : numpy array
        time in hours (e.g., [0., 0.5, 1.0, ...])
    save : bool, optional
        if True, saves the plot as a pdf file, by default False

    Returns
    -------
    plotly figure
    """

    # get index from data and separate it into numbers and letters, save only a unique list of both
    index = data.index.to_list()
    letters = list(set([i[0] for i in index]))
    numbers =list(set([int(i[1:]) for i in index]))

    # sort the lists
    letters.sort()
    numbers.sort()

    # max of y axis
    max_y = data.max(axis=1).max()

    # steps between 0 and max, rounded to 1 decimal
    step = round(max_y/3, 1)
    # this solves a bug when step is too small to be rounded to 1 decimal
    if step == 0.0:
        step = round(max_y/3, 2)

    # max x axis
    max_x = data.columns.max() 
    step_x = round(max_x/2, 0)

    let_len = len(letters)
    num_len = len(numbers)

    fig = make_subplots(rows=let_len, cols=num_len, 
                        shared_xaxes=True, 
                        shared_yaxes=True, 
                        subplot_titles=numbers, 
                        vertical_spacing=0.03, 
                        horizontal_spacing=0.011)

    for i, row in enumerate(letters):
        for j, col in enumerate(numbers):
            # if a combination of row and column is not in the dataframe, skip it
            if f'{row}{col}' not in data.index:
                continue
            else:
                # set color to black and make plots wider
                fig.add_trace(go.Scatter(x=time_h, 
                                        y=data.loc[f'{row}{col}'], 
                                        line=dict(color='black', width=1)), 
                                        row=i+1, col=j+1)
                # y axis update axes, range between min and max of the data
                if j == 0:
                    fig.update_yaxes(range=[0, max_y], gridcolor='white',
                                    title_text=row, 
                                    showline=True, linewidth=1, linecolor='black',mirror=True,
                                    row=i+1, col=j+1, 
                                    tickvals=np.arange(0, max_y, step))
                else:
                    fig.update_yaxes(range=[0, max_y], gridcolor='white',
                                    showline=True, linewidth=1, linecolor='black',mirror=True,
                                    row=i+1, col=j+1, 
                                    tickvals=np.arange(0, max_y, step))
                # x axis, rotate the tick labels by 90 degrees
                fig.update_xaxes(range=[0, 24], row=i+1, col=j+1, gridcolor='white',
                                showline=True, linewidth=1, linecolor='black', mirror=True,
                                tickvals=np.arange(0, max_x+(max_x*0.1), step_x), 
                                ticktext=np.arange(0, max_x+(max_x*0.1), step_x), tickangle=0)
      
           
    # make individual subplots wider 
    # rotate y title_text by 90 degrees
    fig.update_layout(title_text=title, 
                title_x=0.5, title_y=0.95, title_font_size=20,
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                showlegend=False, height=700, width=1000)

    # save the plot in pdf format
    if save:
        fig.write_image(f'{out_name}.pdf')


def plotly_wrapper(time_data, plate, data_type, output):
    """
    Wrapper function to plot the timeseries data using plotly
    Parameters
    ----------
    time_data : pandas dataframe
        The timeseries data to be plotted
    plate : str
        The plate name
    data_type : str
        The data type
    time_h : list
        The time in hours
    Returns
    -------
    None.
    """

    ts = time_data[(time_data.File == plate) & (time_data.Data == data_type)]

    # check if the end of plate is = '.txt' and remove it  
    if plate[-4:] == '.txt':
        plate = plate[:-4] 
    
    # remove non-numeric columns
    time_h = [int(i) for i in time_data.columns if is_number(i)]

    ts_col = time_h.copy()
    ts_col.insert(0, 'Well')

    ts = ts[ts_col]
    ts = ts.set_index('Well')

    time_h = sorted(time_h)
    time_h = np.array(time_h)/60/60

    out_file = f'{output}/{plate}_{data_type}'

    plot_individual_plate_plotly(ts, plate + data_type, out_file, time_h = time_h, save=True)


# filter a list for numeric values
def is_number(s):
    """
    Function to check if a string is a number

    Parameters
    ----------
    s : str
        The string to be checked

    Returns
    -------
    bool
        True if the string is a number, False otherwise
        
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def plot_boxplots(data, grouping_var, temp_var, x_var, y_var, output_dir):
    """
        Function to make boxplots from the dataframe

        Parameters 
        ----------
        data : pandas dataframe
            dataframe containing the data to be plotted
        grouping_var : str
            main variable by which you want to plot the data (e.g. Strains)
        temp_var : str
            variable to be used as a temporary variable to group the data (e.g. the specific strain)
        x_var : str
            variable to be used as the x-axis (e.g. the metformin concentration)
        y_var : str
            variable to be used as the y-axis (e.g. the AUC)
        output_dir : str
            directory to save the plots

        Returns
        -------
        Saves a plot in the specified directory
    """
    
    data = data[data[grouping_var] == temp_var]
    sns.boxplot(x=x_var, y=y_var, data=data)
    sns.swarmplot(x=x_var, y=y_var, data=data, color='black')
    plt.title(temp_var)
    plt.savefig(output_dir + temp_var + '.pdf')
    plt.close()


def get_timeseries_yvar(data, ending):
    """
        Function to get the y_var from a Timeseries dataframe

        Parameters
        ----------
        data : pandas dataframe
            dataframe containing the data to be plotted
        ending : str
            ending of the variable to be used as the y-axis (e.g. 'nm_f' to get the AUC)

        Returns
        -------
        y_var : str
            variable to be used as the y-axis 
    """
    return data[data['Data'].str.endswith(ending)].Data.unique().tolist()[0] 


def plot_simple_timeseries(data, condition, y_var, x_var='Time', title=None):
    """
        Function to plot a timeseries

        Parameters
        ----------
        data : pandas dataframe
            dataframe containing the data to be plotted
        condition : list
            list of variables to be used as the condition
        y_var : str
            variable to be used as the y-axis
        x_var : str
            variable to be used as the x-axis
        title : str
            title of the plot
        
        Returns
        -------
        None

    """
    # set the figure size
    plt.figure(figsize=(8, 7))
    # plot the data, use virids palette
    sns.lineplot(data=data, x=x_var, y=y_var, hue=condition, style=condition, dashes=False, palette='viridis')
    # set the x and y labels
    plt.xlabel('Time (h)')
    plt.ylabel('OD')
    # plot title
    plt.title(title)

# TODO: start separating functions by classes (reading files, converting, analysing...)
