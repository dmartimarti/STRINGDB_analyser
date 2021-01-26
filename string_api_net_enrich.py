#!/usr/bin/env python3

################################################################
## Loads a txt file with genes in a list (in a single column)
## and gets an image file and an enrichment file from STRING
################################################################

from time import sleep
import json
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


def get_net_image(genes,species=511145,output='full_network.svg'):
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

    file_name = output
    print(f"Saving interaction network to {output} file")

    with open(file_name, 'wb') as fh:
        fh.write(response.content)

    sleep(1)


def get_enrichment_data(genes):
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
        "species" : 511145, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org" # your app name

    }

    ## Call STRING
    response = requests.post(request_url, data=params)
    
    # Read the data
    data = json.loads(response.text)
    
    # transform data to a dataframe
    data_long = pd.DataFrame(data)
    
    return data_long


def main():
    genes = gene_list('dnaK_UP_distinct.txt')
    get_net_image(genes)
    enrich = get_enrichment_data(genes)
    sns.countplot(x="category", data=enrich)
    plt.savefig('categories_genes.pdf')

    with pd.ExcelWriter('output.xlsx') as writer:
        for element in enrich.category.unique():
            DF = enrich[enrich['category']==element]
            DF.to_excel(writer, sheet_name=element)


if __name__ == '__main__':
    main()