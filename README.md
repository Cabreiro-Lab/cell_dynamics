# Cellular dynamics 

This script is for internal use for the Cabreiro lab. As we are heavy users of the cell reader Biospa, we are constantly analysing the output of the machine and calculating different parameters of cell growth such as AUCs or the dynamics from the timeseries. 

This script is essentially a rework from the [original script](https://github.com/PNorvaisas/Growth_analysis) used in the lab by [Pov](https://github.com/PNorvaisas), adapted to Python 3.X and with some extra features (multiprocessing, multidimensional analyses, etc). The script is still **under development** and it is not fully tested, any feedback is welcome.

## Installation

The script is written in Python 3.10 and it is recommended to use a virtual environment. The script uses the following libraries:

- numpy
- pandas
- plotly
- tqdm
- itertools
- matplotlib

I highly recommed to create a new conda environment and install the libraries with the following command:

```bash
conda create -n cell_dynamics python=3.10 numpy pandas plotly tqdm matplotlib
```

I will create a requirements.txt file in the future to make the installation easier.

## Usage

The main script is named `growth.py` and depends on the `functions.py` file. The script is called with the following command:

```bash
python growth.py -i <input_folder> -o <output_folder> -t <threads> 
```

The script can be called with the following arguments:

- `-i` or `--input`: Path to the input folder. The script will look for a file named `Design.xlsx` and read it. 
- `-o` or `--output`: Path to the output folder. The script will create a folder for plots and another for the csv files.
- `-t` or `--threads`: Number of threads to use. 


