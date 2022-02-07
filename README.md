# STRINGDB analyser
A quick way to analyse gene/protein sets to investigate networks and functional enrichment.

It comes in two different flavors:
* **string_api_net_enrich.py**: This scripts inputs a list of genes from _E. coli_ or humans in txt format (each gene or protein in a new line), and outputs a network image in svg, a summary of the categories from the different enrichments, and an excel file with all the enrichment information. It will also create radar plots for the most common words for each enrichment category (if any). The file _example.txt_ has a list of genes that can be used to test the script. 
* **string_api_MULTI.py**: This scripts inputs a list of genes from _E. coli_ or humans in excel format. This script is intended to use with gene/protein lists that are UP-regulated and DOWN-regulated. Each direction must be in different Excel sheets, with a sheet name finised in _\_UP'_ or _\_DOWN'_. For example, _sample1\_UP_ and _sample1\_DOWN_. More than one sample can be included in the same Excel file, the script will save each sample in separate subfolders with specific names. The list of genes/proteins **MUST** have a header named _genes_. The file _multi_test_ serves as an example for this script. The script will analyse these samples and directions in a very similar way as the simple script, but creating subfolders for each sample within the Excel file. 



## Requirements

This script runs in Python >3.6 and requires the following libraries: `pandas`, `requests`, `seaborn`,`numpy`, `openpyxl` and `matplotlib`.

There are several ways to install these dependencies. You can either use pip

```bash
pip install pandas
```

Or my favourite way to do these things, through a conda env:

```bash
conda install pandas
```

## Example of use

Go to the folder where you should have both the script and the file you want to analyse (for convenience), and type the following:

```bash
python string_api_net_enrich.py example.txt out_folder ecoli
```

And for the other script, the same:

```bash
python string_api_MULTI.py multi_test.xlsx out_folder ecoli
```

Change input and output for your input file, and your desired output filename.
Right now it allows to specify either _E. coli_ or human as species (type ecoli or human respectively). 

This is an example of the radar plots it's able to extract:

![alt text](https://github.com/dmartimarti/STRINGDB_analyser/blob/main/figs/radar_example.JPG)


## Google colab version

There is a version in Google Colab if you are a bit lazy

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1pMRsXxKZX_L20ehAleMw2yHH79q1uUIK?usp=sharing)


## To do

* _generalise funtions into classes_
* _include more analyses and plots (heatmaps, semantic space of GO terms...)_
* _make word frequency more smart -> concept over words_
* _[long term] build all functions in different files, tidy everything_

### Next action

To adapt the script to two cases (usually up and down regulation) 
