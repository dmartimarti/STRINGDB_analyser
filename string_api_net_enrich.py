#!/usr/bin/env python3

################################################################
## Loads a txt file with genes in a list (in a single column)
## and gets an image file and an enrichment file from STRING
################################################################

from time import sleep
import json
import argparse
import pandas as pd
import requests
import seaborn as sns
import matplotlib.pyplot as plt

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
    string_api_url = "https://string-db.org/api"
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

    with open(f'{out_net}.svg','wb') as fh:
        fh.write(response.content)

    sleep(1)


def get_enrichment_data(genes,species=511145):
    '''
    Function gets gene list and extracts functional enrichment (if any)
    '''
    string_api_url = "https://string-db.org/api"
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

input_file = args.Input
output = args.Output

def main():
    '''
    Main function program
    '''
    genes = gene_list(input_file)
    print(f'Processing file {input_file} with ' + str(len(genes)) + ' elements')
    get_net_image(genes,out_net=output,species=spc)
    enrich = get_enrichment_data(genes,species=spc)
    # plot categories
    fig, ax = plt.subplots()
    g=sns.countplot(x="category", data=enrich)
    g = g.set_xticklabels(g.get_xticklabels() ,rotation=45,
                     horizontalalignment='right')
    fig.tight_layout()
    plt.savefig(f'{output}_categories_enrich.pdf')

    print(f'Saving enrichment in file {output}_output.xlsx')
    with pd.ExcelWriter(f'{output}_output.xlsx') as writer:
        for element in enrich.category.unique():
            enrich_df = enrich[enrich['category']==element]
            enrich_df.to_excel(writer, sheet_name=element)


if __name__ == '__main__':
    main()
