#!/usr/bin/env python3

"""This script analsyses a set of gene/proteins with STRING DB and
gets back a set of plots and tables useful for posterior analyses."""

# Author: Daniel Martinez-Martinez
# Year: 2021 (Pandemic times)

################################################################
## Loads a txt file with genes in a list (in a single column)
## and gets an image file and an enrichment file from STRING
################################################################

from time import sleep
from collections import Counter

import os
import json
import argparse
import pandas as pd
import requests
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# define functions

def gene_list(file_list):
    '''
    reads a gene list from a txt files
    '''
    data = pd.read_csv(file_list, sep=' ', header=None)
    genes = data[0].values.tolist()
    return genes


def get_net_image(genes,species=511145,out_net='full_network.svg'):
    '''
    This function gets a gene list as an input and
    outputs a svg image of the network from those genes
    from strin- db
    A different species and output name can be chosen
    '''
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "svg"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    ## Parameters
    params = {
            "identifiers" : "\r".join(genes), # your protein
            "species" : species, # species NCBI identifier
            "network_flavor": "confidence", # show confidence links
    }

    response = requests.post(request_url, data=params)

    print(f"Saving interaction network to {out_net}.svg file")

    with open(f'./{output}/{out_net}.svg','wb') as fh_net:
        fh_net.write(response.content)

    sleep(1)


def get_enrichment_data(genes,species=511145):
    '''
    Function gets gene list and extracts functional enrichment (if any)
    '''
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"

    ## Construct the request
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    params = {

        "identifiers" : "%0d".join(genes), # your protein
        "species" : species, # species NCBI identifier
        "caller_identity" : "www.awesome_app.org" # your app name

    }

    ## Call STRING
    response = requests.post(request_url, data=params)
    # Read the data
    data = json.loads(response.text)
    # transform data to a dataframe
    data_long = pd.DataFrame(data)
    return data_long


def count_words(df,category='Process',nwords=10):
    '''
    This function inputs the enrichment dataframe from get_enrichment_data
    and outputs a list of names with their relative frequencies
    It can be parsed per category type
    The function controls for useless words (e.g. to, and, of...) and for 
    other punctuation symbols
    '''
    proc_names = df[df['category']==category]['description'].tolist()
    words = []
    for elm in proc_names:
        for word in elm.split():
            word = word.replace(',', '')
            words.append(word.lower())
    # create a data frame
    word_df = pd.DataFrame.from_dict(Counter(words), orient='index',columns=['count'])
    # filter words and keep the top ones
    filter_list = ['process', 'substance', 'to', 'a','metabolic','via',
                    'and','of','incl.','by','in','with','the','from']
    word_df = word_df[~word_df.index.isin(filter_list)]
    word_df = word_df.sort_values('count',ascending=False)
    # calculate the relative word use in this sublist
    word_df = word_df.head(nwords)
    word_df = word_df.assign(relative= word_df['count']/ word_df['count'].sum())
    return word_df


def radar_chart(dft,category):
    '''
    This function inputs the word count resulting from count_words
    and plots the radar chart from that list

    '''

    # transpose the dataframe
    dft = dft.transpose()

    # Each attribute we'll plot in the radar chart.
    labels = list(dft)[0:]
    # Let's take the values
    values = dft.loc['relative'].tolist()

    # Number of variables we're plotting.
    num_vars = len(labels)

    # Split the circle into even parts and save the angles
    # so we know where to put each axis.
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # The plot is a circle, so we need to "complete the loop"
    # and append the start value to the end.
    values += values[:1]
    angles += angles[:1]

    # ax = plt.subplot(polar=True)
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    # Draw the outline of our data.
    ax.plot(angles, values, color='red', linewidth=1)
    # Fill it in.
    ax.fill(angles, values, color='red', alpha=0.25)

    # Fix axis to go in the right order and start at 12 o'clock.
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Draw axis lines for each angle and label.
    ax.set_thetagrids(np.degrees(angles[:-1]), labels)

    # Go through labels and adjust alignment based on where
    # it is in the circle.
    for label, angle in zip(ax.get_xticklabels(), angles):
        if angle in (0, np.pi):
            label.set_horizontalalignment('center')
        elif 0 < angle < np.pi:
            label.set_horizontalalignment('left')
        else:
            label.set_horizontalalignment('right')

    ax.set_rlabel_position(180 / num_vars) # sets y-label to middle

    # Add some custom styling.
    ax.tick_params(colors='#222222') # Change the color of the tick labels.
    ax.tick_params(axis='y', labelsize=8) # Make the y-axis (0-100) labels smaller.
    ax.grid(color='#AAAAAA') # Change the color of the circular gridlines.
    ax.spines['polar'].set_color('#222222') # Change color of the outermost gridline
    ax.set_facecolor('#FAFAFA') # circle background color
    fig.tight_layout()
    # Lastly, give the chart a title and give it some
    # padding above the "Acceleration" label.
    ax.set_title(category, y=1.08)

# define options to parse
# Create the parser
my_parser = argparse.ArgumentParser(
    prog='STRING API enrich',
    description='List the content of a folder')

# Add the arguments
my_parser.add_argument('Input',
                       metavar='-i',
                       type=str,
                       help='input file')
my_parser.add_argument('Output',
                       metavar='-o',
                       type=str,
                       help='output name for output files')
my_parser.add_argument('Species',
                       metavar='-s',
                       type=str,
                       help='select between ecoli or human')

# Execute the parse_args() method
args = my_parser.parse_args()

# redefine args
input_file = args.Input
output = args.Output

# try creation of folder
try:
    os.mkdir(output)
except OSError:
    print (f"Creation of the directory {output} failed")
else:
    print (f"Successfully created the directory {output}")

# program version
_VERSION_ = 0.3

# define the list of species included in the script
species_list = {
    'ecoli':511145,
    'human':9606
}

# define species for the analysis
spc = species_list[args.Species]


def main():
    '''
    Main function program
    '''
    genes = gene_list(input_file)
    print(f'Processing file {input_file} with ' + str(len(genes)) + ' elements')
    get_net_image(genes,out_net=output,species=spc)
    enrich = get_enrichment_data(genes,species=spc)
    # plot categories

    # check that the enrich is not emtpy
    if not enrich.empty:
        if len(enrich.category.unique()) > 0:
            # print summary of categories
            fig, ax = plt.subplots()
            g_plot=sns.countplot(x="category", data=enrich)
            g_plot = g_plot.set_xticklabels(g_plot.get_xticklabels() ,rotation=45,
                             horizontalalignment='right')
            fig.tight_layout()
            plt.savefig(f'./{output}/{output}_categories_enrich.pdf')

            print(f'Printing radar plots for {enrich.category.unique()}')
            for cat in enrich.category.unique():
                word_df = count_words(enrich,category=cat,nwords=10)
                radar_chart(word_df,category=cat)
                plt.savefig(f'./{output}/{cat}_radar_chart.pdf')
        else:
            print('There are not categories to plot')


        print(f'Saving enrichment in file {output}_output.xlsx')
        with pd.ExcelWriter(f'./{output}/{output}_output.xlsx') as writer:
            for element in enrich.category.unique():
                enrich_df = enrich[enrich['category']==element]
                enrich_df.to_excel(writer, sheet_name=element)

    elif enrich.empty:
        print('There were no enriched categories!')

    print(f'All analyses have been finished for the file {input_file}')
    
if __name__ == '__main__':
    main()