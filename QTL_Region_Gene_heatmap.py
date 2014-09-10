#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
from urllib2 import urlopen
import numpy as np
import pylab as pl
import matplotlib.colors as colors
import random
from matplotlib import gridspec

def usage():
    test="name"
    message='''
python QTL_Gene_heatmap.py --input chr3.headingdate.ggplot.table --rnaseq chr3.headingdate.ggplot.ColdStress.table  --output chr3.headingdate.ggplot

    '''
    print message

def set_ticks_X(ax, ylen, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 8)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(qtl.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(ylab, minor=False, fontsize=8)
    ax.set_yticklabels([])

    # rotate the
    plt.xticks(rotation=40)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def set_ticks_LOD_X(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 8)

    # turn off the frame
    #ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(qtl.shape[0]) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    #ax.set_xticklabels([])
    ax.set_yticklabels([])

    # rotate the
    plt.xticks(rotation=40, fontsize=8)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    #for t in ax.xaxis.get_major_ticks():
    #    t.tick1On = False
    #    t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def set_ticks_XY(ax, xlen, ylen, xlab, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(xlab, minor=False)
    ax.set_yticklabels(ylab, minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def set_ticks_XY_Right(ax, xlen, ylen, xlab, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(xlab, minor=False, fontsize=8)
    ax.set_yticklabels(ylab, minor=False)

    # Set the y label position relative to y axis 
    ax.tick_params(axis='y',direction='out', pad=140)

    # rotate the
    plt.xticks(rotation=40)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def linecount(filename):
    return len(open(filename).readlines())

def colorbar(heatmap,  cbaxes, label, ticks, ticke):
    # Add a colorbar below the heatmap
    # add_axes refer to [left, bottom, width, height], where the coordinates are just fractions that go from 0 to 1 of the plotting area.
    #cbaxes = fig.add_axes([0.15, 0.15, 0.2, 0.02]) 
    cb = plt.colorbar(heatmap, cax = cbaxes, orientation='horizontal', ticks=np.linspace(ticks, ticke, 5)) 
    #cb = pl.colorbar(heatmap, orientation='horizontal', shrink=0.5825,
    #                 fraction=0.02, pad=-0.035,ticks=np.linspace(-1, 1, 5),
    #                 use_gridspec=True)
    cb.set_label('$\mathrm{%s}\ $' %(label))
    #return fig

def heatmap(gene_list, rnaseq1_list, rnaseq2_list, rnaseq3_list, annotation, lod, project):
    
    genecount = int(linecount(gene_list)) - 1 
    fheight = genecount * 0.25
    titley  = 1
    colorb1  = []
    colorb2  = []
    if genecount < 90:
        titley += 0.06
        colorb1 = [0.12, 0.17, 0.2, 0.01]
        colorb2 = [0.47, 0.17, 0.2, 0.01]            
    elif genecount < 200:
        titley += 0.03
        colorb1 = [0.12, 0.18, 0.2, 0.005]
        colorb2 = [0.47, 0.18, 0.2, 0.005]
    else:
        titley += 0.01
        colorb1 = [0.12, 0.195, 0.2, 0.0015]
        colorb2 = [0.47, 0.195, 0.2, 0.0015]

    # Create a figure.
    figsize=(8,fheight)
    fig = pl.figure(figsize=figsize)
    pl.subplots_adjust(bottom=0.2)
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.5, 3, 3, 3, 3])

    # Read Dataframe from file
    qtl = pd.read_table(gene_list, index_col=0)
    rnaseq1 = pd.read_table(rnaseq1_list, index_col=0)
    rnaseq2 = pd.read_table(rnaseq2_list, index_col=0)
    rnaseq3 = pd.read_table(rnaseq3_list, index_col=0)
    anno = pd.read_table(annotation, sep='\t')
    lod  = pd.read_table(lod, index_col=0)

    #use truncated color table from 0.2 to 0.8
    cmap = plt.cm.get_cmap("Blues")
    new_cmap = truncate_colormap(cmap, 0.2, 0.8)
 

    # Draw heatmap for variations
    #ax = fig.add_subplot(131, frame_on=False)
    #ax = fig.add_subplot(133, frame_on=True)
    ax0 = plt.subplot(gs[0])
    ax0.set_ylim(0, genecount)
    heatmap = ax0.pcolor(qtl, cmap=new_cmap, alpha=0.8, vmin= 0, vmax= 4)
    # Set names for x,y_ticks
    ax0 = set_ticks_XY(ax0, qtl.shape[0], qtl.shape[1], qtl.columns, qtl.index)
    
    ticks0 = 0
    ticke0 = 4 
    #cbaxes0 = fig.add_axes([0.12, 0.18, 0.2, 0.01])
    cbaxes0 = fig.add_axes(colorb1)
    colorbar(heatmap, cbaxes0, 'Variation Effect', ticks0, ticke0)
    
    # Draw heatmap for gene expression: cold
    #ax = fig.add_subplot(132, frame_on=False) 
    ax1 = plt.subplot(gs[1])
    ax1.set_ylim(0,genecount)
    heatmap = ax1.pcolor(rnaseq1, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)  
    plt.title('Cold', y= titley) 
    # Set names for x_ticks only
    ax1 = set_ticks_X(ax1, rnaseq1.shape[1], rnaseq1.columns)
    #ax = fig.add_subplot(133, frame_on=True)
    
    ticks1 = -2
    ticke1 = 2
    #cbaxes1 = fig.add_axes([0.47, 0.18, 0.2, 0.01])
    cbaxes1 = fig.add_axes(colorb2)
    colorbar(heatmap, cbaxes1, 'Log2(HEG4/NB)', ticks1, ticke1)

   # Draw heatmap for gene expression: drought
    #ax = fig.add_subplot(132, frame_on=False)
    ax2 = plt.subplot(gs[2])
    ax2.set_ylim(0, genecount)
    heatmap = ax2.pcolor(rnaseq2, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)
    plt.title('Drought', y = titley) 
    # Set names for x_ticks only
    ax2 = set_ticks_X(ax2, rnaseq2.shape[1],  rnaseq2.columns)

    # Draw heatmap for gene expression: salt
    #ax = fig.add_subplot(132, frame_on=False)
    ax3 = plt.subplot(gs[3])
    ax3.set_ylim(0, genecount)
    heatmap = ax3.pcolor(rnaseq3, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)
    plt.title('Salt', y= titley) 
    # Set names for x_ticks only
    #ax3 = set_ticks_X(ax3, rnaseq3.shape[1],  rnaseq3.columns)
    ylabs = anno['Annotation']
    ax3 = set_ticks_XY_Right(ax3, rnaseq3.shape[0], rnaseq3.shape[1], rnaseq3.columns, ylabs)
    
 
    # Draw LOD curve
    ax4 = plt.subplot(gs[4])
    x=lod['LOD']
    y=np.arange(lod.shape[0]) + 0.5
    ax4.set_xlim(min(x)*0.8, max(x)*1.2)
    ax4.set_ylim(0, genecount)
    plt.title('LOD', y= titley)
    ax4.plot(x,y)
    # Set names for x_ticks 
    ax4 = set_ticks_LOD_X(ax4)
    #ylabs = anno['Annotation']
    #ax3 = set_ticks_XY_Right(ax3, rnaseq3.shape[0], rnaseq3.shape[1], rnaseq3.columns, ylabs)

    # Save file
    fig.savefig('%s.pdf' %(project), bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--rnaseq1')
    parser.add_argument('--rnaseq2')
    parser.add_argument('--rnaseq3')
    parser.add_argument('--anno')
    parser.add_argument('--lod')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.rnaseq1) > 0
    except:
        usage()
        sys.exit(2)

    heatmap(args.input, args.rnaseq1, args.rnaseq2, args.rnaseq3, args.anno, args.lod, args.output)

if __name__ == '__main__':
    main()

