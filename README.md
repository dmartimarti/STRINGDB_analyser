# STRINGDB analyser
A quick way to analyse gene/protein sets to investigate networks and functional enrichment.

This scripts inputs a list of genes from E. coli in txt format (each gene or protein in a new line), and outputs a network image in svg, a summary of the categories from the different enrichments, and an excel file with all the enrichment information. 

## Requirements

This script runs in Python >3.6 and requires the following libraries: `pandas`, `requests`, `seaborn` and `matplotlib`.

If you are a Mac user, go to the Terminal and type:

```bash
pip install pandas
```

If you are a Windows user, I recommed to install Anaconda as your environment manager, go to the taskbar and search for `Anaconda prompt`. Once there, just install the libraries as follows:

```bash
pip install pandas
```

or

```bash
conda install pandas
```

## Example of use

Go to the folder where you should have both the script and the file you want to analyse (for convenience), and type the following:

```bash
python string_api_net_enrich.py [input] [output] [species]
```

Change input and output for your input file, and your desired output filename.
Right now it allows to specify either _E. coli_ or human as species (type ecoli or human respectively). 
