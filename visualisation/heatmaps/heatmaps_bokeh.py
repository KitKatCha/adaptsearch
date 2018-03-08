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
import argparse
from math import pi
import logging
logging.basicConfig()

def makeLineOrderFromTree(tree):

    t = open(tree, "r")
    arbre = t.readline()
    t.close()

    mapping = [ ('(', ''), (')', ''), (':', ''), (';', ''), ('.', ''), ('#', '')]
    for k, v in mapping:
        arbre = arbre.replace(k, v)

    arbre = arbre.split(',')
    arbre = [sp[0:2] for sp in arbre]

    return arbre


def makeHeatmap(df, what, on_what, colors, tree, range_x):

    # Select the good lines
    df = {
        'countings': df[0::3],
        'frequencies': df[1::3],
        'pvalues': df[2::3]
    }[what]

    y = list(df.index)
    x = list(df.columns)
    
    # heatmap width - apparently, width cannot be smaller than height/2
    w = {
        4: 200,
        20: 400,
        61: 900
    }[len(x)]

    # heatmap height
    # TODO

    # Keep only species name abbreviation
    name_index = {elem:str.split(elem, "_")[0] for elem in y}
    df = df.rename(index=name_index)
    df = pd.DataFrame(df.stack(), columns=['counts']).reset_index()

    # Lines order according to phylogenetic tree    
    range_names = makeLineOrderFromTree(tree)
    # title
    title = "{} of {}".format(what, on_what)
    # Color range    
    mapper = LinearColorMapper(palette=colors, low=df.counts.min(), high=df.counts.max())
    
    # format to plot
    print df
    source = ColumnDataSource(df)
    # Heatmap
    p = figure(title=title,
           x_range=range_x, y_range=range_names,
           x_axis_location="above", plot_width=w, plot_height=400,
           toolbar_location=None, tools=""
           )

    # Misc. styling
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "10pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3

    # Draw rectangles
    p.rect(x=on_what, y="Species", width=1, height=1,
           source=source,
           fill_color={'field': 'counts', 'transform': mapper},
           line_color=None)
    """
    p.circle(x=on_what, y="Species", size=15,
       source=source,
       fill_color={'field': 'counts', 'transform': mapper},
       line_color=None)
    """

    # Add color bar
    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="10pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format="%.3f"),
                         label_standoff=15, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')

    export_png(p, filename='{}.png'.format(title.replace(" ", "_")))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input csv file")
    parser.add_argument("represent", choices=("countings", "frequencies", "pvalues", "all"), help="choose what to plot")
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
    
    pvals = df[2::3]
    y = list(df.index)
    name_index = {elem:str.split(elem, "_")[0] for elem in y}
    pvals = pvals.rename(index=name_index)
    pvals = pd.DataFrame(pvals.stack(), columns=['pvalues']).reset_index()
    print pvals

    if args.represent != "all":
        makeHeatmap(df, args.represent, on_what, colors, args.tree, x_range[on_what])
    elif args.represent == "all":
        what = ("countings", "frequencies", "pvalues")
        makeHeatmap(df, what[0], on_what, colors, args.tree, x_range[on_what])
        makeHeatmap(df, what[1], on_what, colors, args.tree, x_range[on_what])
        makeHeatmap(df, what[2], on_what, colors, args.tree, x_range[on_what])

    # TODO :
    #   - taille des rectangles et de la figure en fonction du nb de colonnes
    #   - sort des colonnes : alphabétique / par groupe
    #   - sort des lignes : alphabétique / groupe
    
if __name__ == "__main__":
    main()