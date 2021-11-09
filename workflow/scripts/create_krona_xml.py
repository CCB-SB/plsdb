#!/usr/bin/python

"""
Krona XML 2.0 specification: https://github.com/marbl/Krona/wiki/Krona-2.0-XML-Specification
<krona ...>
    <attributes magnitude="sample_attr_1">
        <attribute ...>sample_attr_1</attribute>
        <attribute ...>sample_attr_2</attribute>
        <list>sample_list_1</list>
    </attributes>
    <color ...></color>
    <datasets>
        <dataset>set1</dataset>
        <dataset>set2</dataset>
    </datasets>
    <node name="root">
        <sample_attr_1>
            <val>7</val>
            <val>5</val>
        </sample_attr_1>
        <sample_attr_2>...</sample_attr_2>
        <sample_list_1>
            <vals>
                <val>12</val>
                <val>15</val>
            </vals>
            <vals>...</vals>
        </sample_list_1>
        <node name="child">...</node>
    </node>
</krona>
"""

import pandas
import logging
import argparse
import xml.etree.cElementTree as ET

from utils import setup_logger


##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tab', '-t', help="Plasmid info table", required=True)
    parser.add_argument('--ofile', '-o', help="Output file (*.xml)", required=True)
    return parser

##################################################
# FUNC
##################################################
def new_taxon_name(taxon):
    t_name, t_id = taxon
    if pandas.isnull(t_name):
        assert pandas.isnull(t_id)
        return 'unknown'
    return '{} (ID {})'.format(t_name, int(t_id))

def aggr_taxa(node, df, rank):
    # print('Rank %s, %d entries' % (rank, df.shape[0]))
    for taxon, taxon_df in df.groupby(by='new_%s' % rank):
        taxon_node = ET.SubElement(node, 'node')
        taxon_node.attrib['name'] = taxon
        taxon_count = ET.SubElement(taxon_node, 'count')
        for db in labels:
            if db == 'All':
                ET.SubElement(taxon_count, 'val').text = str(taxon_df.shape[0])
            else:
                ET.SubElement(taxon_count, 'val').text = str(taxon_df.loc[taxon_df['Source_NUCCORE'] == db,:].shape[0])
        if ranks.index(rank) < len(ranks) - 1:
            taxon_node = aggr_taxa(taxon_node, taxon_df, ranks[ranks.index(rank) + 1])
    return node

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    args = get_arg_parser().parse_args()

    # taxonomy ranks
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # read in data
    tab = pandas.read_csv(args.tab, sep='\t', header=0)

    # labels
    labels = ['All'] + sorted(set(tab['Source_NUCCORE']))

    # create new taxon names: "taxon <name> (<taxon ID>)"
    for rank in ranks:
        tab['new_%s' % rank] = tab[['taxon_%s_name' % rank, 'taxon_%s_id' % rank]].apply(new_taxon_name, axis=1)
        logging.info('Processed {} taxa'.format(rank))

    ##############################
    # XML START
    logging.info('Create XML tree structure...')

    # root element
    krona = ET.Element("krona")

    # attributes
    attrs = ET.SubElement(krona, 'attributes')
    attrs.attrib['magnitude'] = "count"
    attr_count = ET.SubElement(attrs, 'attribute')
    attr_count.attrib['display'] = "Count"
    attr_count.text = "count"

    # datasets
    datasets = ET.SubElement(krona, 'datasets')
    for db in labels:
        ET.SubElement(datasets, 'dataset').text = db

    # taxonomy (nodes)
    # root node -> all
    root_node = ET.SubElement(krona, 'node')
    root_node.attrib['name'] = "cellular organisms"
    root_count = ET.SubElement(root_node, 'count')
    for db in labels:
        if db == 'All':
            ET.SubElement(root_count, 'val').text = str(tab.shape[0])
        else:
            ET.SubElement(root_count, 'val').text = str(tab.loc[tab['Source_NUCCORE'] == db,:].shape[0])
    # nodes for ranks
    nodes = aggr_taxa(root_node, tab, 'superkingdom')

    # XML END
    ##############################

    # create tree
    tree = ET.ElementTree(krona)
    # save to file
    tree.write(args.ofile)
