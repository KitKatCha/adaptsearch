#!/usr/bin/env python
# coding: utf8

from bokeh.io import show, output_file, export_png

from bokeh.models import (
    ColumnDataSource,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar
)
from bokeh.plotting import figure

import pandas as pd
import numpy as np
import argparse
from math import pi
import logging
logging.basicConfig()

def makeLineOrderFromTree(tree):
    """ Read a newick tree and return a list of the species, ordered as they come in the tree

    Args : tree (file) : a newick tree
    Return (list of str)

    """

    t = open(tree, "r")
    arbre = t.readline()
    t.close()

    mapping = [ ('(', ''), (')', ''), (':', ''), (';', ''), ('.', ''), ('#', '')]
    for k, v in mapping:
        arbre = arbre.replace(k, v)

    arbre = arbre.split(',')
    arbre = [sp[0:2] for sp in arbre]

    return arbre

def makeHeatmap(df, what, on_what, colors, range_x, range_y, width, height):
    """ Make a bokeh heatmap from a pandas dataframe 

    Args : 
        - df (pd df) : an already stacked pandas dataframe
        - what (str) : indicates the column name being plotted : 'countings' or 'frequencies'
        - on_what (str) : indicates the type of data : 'codons', 'amino-acids', 'amino-acids-types'
        - colors (list of str) : list with hexadecimals colors
        - range_x (list of str) : list of 'on_what' for an ordered x_axis
        - range_y (list of str) : list of the species name for an ordered y_axis. Obtained by makeLineOrderFromTree()
        - width (int) : the figure height
        - height (int) : the figure width

    """
    title = "{} of {}".format(what, on_what)

    # Color range
    # df.what did not work (what was not recognized as a variable)
    mapper = LinearColorMapper(palette=colors, low=df.loc[:,[what]].min()[-1], high=df.loc[:,[what]].max()[-1])
    
    # Heatmap. no tool bar for png
    p = figure(title=title,
           x_range=range_x, y_range=range_y,
           x_axis_location="above", plot_width=width, plot_height=height,
           toolbar_location=None, tools=""
           )

    # Misc. styling
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "10pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3
    
    p.circle(x=on_what, y="Species", size='pvalues',
       source=df,
       fill_color={'field': what, 'transform': mapper},
       line_color=None)

    # legend numbers format
    if what == "countings":
        frmt = "%d"
    elif what == "frequencies":
        frmt = "%.3f"

    # Add color bar
    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="10pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format=frmt),
                         label_standoff=15, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')

    # output
    export_png(p, filename='{}.png'.format(title.replace(" ", "_")))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input csv file")
    parser.add_argument("represent", choices=("countings", "frequencies", "all"), help="choose what to plot")
    parser.add_argument("tree", help="RAxML tree for lines order")
    args = parser.parse_args()

    print "\n*** Heatmaps on the results of MutCount in concatenated mode ***"
    print "\n    Species are ordered according to their position in the tree."
    print "    Codons are ordered according to their corresponding amino-acids,"
    print "    and amino-acids are ordered accordint to their classification."

    df = pd.read_csv(args.file, sep=",", index_col=0)

    # Prepare data    
    on_what = {
        4: "amino-acids-types",
        20: "amino-acids",
        61: "codons"
    }[len(list(df.columns))]

    df.columns.name = on_what

    x = list(df.columns)
    
    # heatmap width - apparently, width cannot be smaller than height/2
    width = {
        4: 300,
        20: 500,
        61: 1000
    }[len(x)]

    y = list(df.index)

    # replacing pvalues by corresponding rect size # Find a one-liner to do that
    # Maybe consider to add another line instead of replacing
    for i in range(2,len(y),3):
        l=df.iloc[i]
        for j in range(0,len(x)):
            if l[j] < 0.025:
                l[j] = 5
            elif l[j] > 0.975:
                l[j] = 10
            else:
                l[j] = 15

    # correct names and stacking
    name_index = {elem:str.split(elem, "_")[0] for elem in y}
    df = df.rename(index=name_index)
    c = np.array(['countings','frequencies','pvalues'])
    df.index = [df.index, c[np.arange(len(df.index)) % 3]]
    df = df.stack().unstack(1).reset_index()

    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]

    # choose order of lines and columns :
    # x_range et y_range dans figure
    x_range = {'codons' : ['ttt','ttc','tgg','tat','tac','gat','gac','gaa','gag','cat',
                            'cac','aaa','aag','cgt','cgc','cga','cgg','aga','agg','tgt',
                            'tgc','aat','aac','cct','cca','ccg','ccc','caa','cag','tct',
                            'tcc','tca','tcg','agt','agc','act','acc','aca','acg','gct',
                            'gcc','gca','gcg','ggt','ggc','gga','ggg','att','atc','ata',
                            'tta','ttg','ctt','ctc','cta','ctg','atg','gtt','gtc','gta','gtg'],
                'amino-acids' : ['F','W','Y','D','E','H','K','R','C','N','P','Q','S','T','A','G','I','L','M','V'],
                'amino-acids-types' : ['aromatics','charged','polar','unpolar']}

    # Lines order according to phylogenetic tree
    y_range = makeLineOrderFromTree(args.tree)
    height = len(y_range)*50

    if args.represent == "countings":
        makeHeatmap(df, args.represent, on_what, colors, x_range[on_what], y_range, width, height)
    elif args.represent == "frequencies":
        makeHeatmap(df, args.represent, on_what, colors, x_range[on_what], y_range, width, height)
    elif args.represent == "all":        
        makeHeatmap(df, "countings", on_what, colors, x_range[on_what], y_range, width, height)
        makeHeatmap(df, "frequencies", on_what, colors, x_range[on_what], y_range, width, height)

    # TODO :   
    #   - cluster a posteriori sur similitude
    #   - hover : valeur comptages/frequence + pvalue -> implique une colonne pvalue en plus de la colonne size
    
if __name__ == "__main__":
    main()