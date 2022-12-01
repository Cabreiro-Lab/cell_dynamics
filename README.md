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
- seaborn
- scipy
- openpyxl

I highly recommed to create a new conda environment and install the libraries with the following command:

```bash
conda create -n <env_name> python=3.10 numpy pandas plotly tqdm matplotlib seaborn scipy openpyxl
```

Activate the new environment.

```bash
conda activate <env_name>
```

As an extra step to save plots in Plotly, you need to install kaleido.

```bash
conda install -c conda-forge python-kaleido
```

*I will create a requirements.txt file in the future to make the installation easier.*

After you have created a new environment and installed all the libraries, **you need to clone the repository in your local machine**.

```bash
cd <path_to_your_local_folder>
git clone https://github.com/Cabreiro-Lab/cell_dynamics
```

The scripts are in the `scripts` folder, which contains the main script (`growth.py`), and a folder with necessary functions. Right now there's only one extra script in there, but there might be more in the future. If you want to run the script, you need to copy `growth.py` and the functions folder to the folder where you have your data. **However, if you want to run it anywhere**, keep reading.

## OPTIONAL: make the script executable anywhere

If you want to make this script available from anywhere, you can add the path to the script to the `PATH` variable in your `.zshrc` file. For example, if you have the github repository in the folder `~/Documents/`, you can add the following line to your `.zshrc` file:

**Be aware that modifying zshrc can produce unexpected results**

```bash
echo "export PATH=$PATH:~/Documents/cell_dynamics/scripts" >> ~/.zshrc
```

Add also the functions folder:
    
```bash
echo "export PATH=$PATH:~/Documents/cell_dynamics/scripts/functions" >> ~/.zshrc
```

Then go to the scripts folder and make the main script executable:

```bash
cd ~/Documents/cell_dynamics/scripts
chmod +x growth.py
```

**If your folders are different, change them accordingly!**

This should make the script available from anywhere in your computer. You might need to restart your terminal to make it work by typing `source ~/.zshrc`, or open a new terminal.

## Usage

The main script is named `growth.py` and depends on the `functions.py` file. The script is called with the following command:

```bash
python growth.py -i <input_folder> -o <output_folder> -t <threads> 
```

The script can be called with the following arguments:

- `-i` or `--input`: Path to the input folder **where Design.xlsx is**. The script will look for a file named `Design.xlsx` and read it. 
- `-o` or `--output`: Path to the output folder. The script will create a folder for plots and another for the csv files.
- `-t` or `--threads`: Number of threads to use. 

### Design.xlsx

A word about the *input file*: for now it needs to be named exactly as `Design.xlsx`, have a sheet named `Design` and it needs to have the following columns:
- `File`: Name of the txt files to analyse. The script will look for a file with the same name in the input folder.

Any other sheet in the input file will be ignored.

If you want to include information about the plate pattern, `Design.xlsx` must have a column named `Pattern`, where it indicates the name of the Pattern file or files that it will read and parse. This pattern file can have as many sheets as you want, with a shape of a 96-well plate. If you don't want a specific column to be read by the script, you can name it starting with an underscore, e.g., `_Variable`.

**New function**, now is possible to include a sheet named `analysis` in the Design.xlsx file. This sheet will be used to specify **two variables**: a *grouping variable* and a *condition*. The grouping variable will be used to group the data, and the condition will be used as a x-axis variable for the boxplot. The script will look for the following columns: `grouping_variable` and `condition`. Below the column names you must specify the name of the column in the `Design` sheet that you want to use as grouping variable and condition. The variables to be used can be variables used in the Design file and/or in the Pattern files. 

For example, imagine you are analysing different strains with 0 and 50 mM of metformin. You can include as a grouping variable `strain`, and as a condition `metformin`. The script will group the data by strain and plot the boxplot for each strain, with the metformin as the x-axis variable.

*TODO*: add more flexibility to the grouping variable, so it can be a combination of variables.

## Output

The script will create a folder named as the specified output, and within it will create a folder named `Plots`. It will save two .csv files in the output folder, one with the AUCs and another with the timeseries. Then it will save all the plots within the `Plots` folder.

It will output 3 types of plots:
- A timeseries plot with raw values that will be named just as the input file name and the wave length, e.g., *plate_1_595nm.pdf*
- A timeseries plot a Wiener smoothing filter applied, named as the input file name and the wave length + '_f', e.g., *plate_1_595nm_f.pdf*
- A plot of the growth rates per time, named as the input file name and the wave length + '_dt', e.g., *plate_1_595nm_dt.pdf*

The same variables can be found in the AUCs and timeseries .csv files as column names. 

## License

MIT License

Copyright (c) [2022] [Daniel Martinez Martinez]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
