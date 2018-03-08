#!/usr/bin/env python
# coding: utf8

from bokeh.io import show, output_file
from bokeh.io import export_png

from bokeh.models import ColumnDataSource, FactorRange
from bokeh.models import LinearAxis, Range1d
from bokeh.models import Legend

from bokeh.palettes import Category20c, magma
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.core.properties import value

import pandas as pd
import argparse
import logging
logging.basicConfig()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input csv file")
    args = parser.parse_args()

    df = pd.read_csv(args.file, sep=",", index_col=0)

    # 1 - species for xaxis 
    list_sp = [str.split(sps, "_")[0] for sps in list(df.index.values)[0::3]]

    # 2 - codons/aa/aatypes names
    subx = list(df)
    """
    if len(subx) == 4:
        palette = ["#c9d9d3", "#718dbf", "#e84d60", "#666dca"]
        end_title = "amino-acids-types"
    if len(subx) == 20:
        palette = Category20c[20]
        end_title = "amino-acids"
    if len(subx) == 61:
        palette = magma(61)
        end_title = "codons"
    """
    palette = {
        4: ["#c9d9d3", "#718dbf", "#e84d60", "#666dca"],
        20: Category20c[20],
        61: magma(61)
    }[len(subx)]
    end_title = {
        4: "amino-acids-types",
        20: "amino-acids",
        61: "codons"
    }[len(subx)]

    # 3 - store values
    data_counts = {}
    data_counts['xaxis'] = list_sp
    for i in range(len(list(df))):        
        data_counts[subx[i]]=list(df[0::3][subx[i]])

    data_freqs = {}
    data_freqs['xaxis'] = list_sp
    for i in range(len(list(df))):        
        data_freqs[subx[i]]=list(df[1::3][subx[i]])

    # 4 - Plot Counts
    title = 'Countings on all species ({})'.format(end_title)    
    source = ColumnDataSource(data=data_counts)
    p = figure(x_range=list_sp, plot_height=500, title=title, toolbar_location=None, tools="")
    p.vbar_stack(subx, x='xaxis', width=0.9, source=source, color=palette, legend=[value(x) for x in subx])
    p.legend.location = "top_right"
    p.legend.orientation = "horizontal"

    # remove legend from inside and placed outside
    new_legend = p.legend[0]
    p.legend[0].plot = None
    p.add_layout(new_legend, 'below')

    export_png(p, filename='{}.png'.format(title.replace(" ", "_")))

    # 5 - Plot freqs
    title = 'Frequencies on all species ({})'.format(end_title)    
    source = ColumnDataSource(data=data_freqs)
    p = figure(x_range=list_sp, plot_height=500, title=title, toolbar_location=None, tools="")
    p.vbar_stack(subx, x='xaxis', width=0.9, source=source, color=palette, legend=[value(x) for x in subx])
    p.legend.location = "top_right"
    p.legend.orientation = "horizontal"

    # remove legend from inside and placed outside
    new_legend = p.legend[0]
    p.legend[0].plot = None
    p.add_layout(new_legend, 'below')

    export_png(p, filename='{}.png'.format(title.replace(" ", "_")))    

if __name__ == "__main__":
    main()