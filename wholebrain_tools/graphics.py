import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from matplotlib.ticker import FixedLocator, MaxNLocator
import seaborn as sns
import pandas as pd
import typing as tp
import copy

from scipy.stats import spearmanr, pearsonr
from statsmodels.stats import anova
from sklearn.linear_model import RANSACRegressor, HuberRegressor

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

def midOntologyHeatmap(data:pd.DataFrame, atlas:object, vmin:float=0, vmax:float=3, 
    cmap:str='PuBu', figsize:tuple=(5,12), fontScaling:float=1, title:str=None):

    """
    Creates a vertical heatmap of regions at mid-ontology level from a dataframe.
    You can customize parts of the plot with some of the additional arguments or
    you can use the returned graphical object to customize it

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        vmin:float=0
            Minimum value for the heatmap
        vmax:float=3
            Maximum value for the heatmap
        cmap:str='PuBu',
            Colormap
        figsize:tuple=(5,12)
            Size of the figure
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure

        Returns
        -------
        axObj:obj
            Graphical object normally returned by seaborn.
    """

    # Set custom color for NaN values
    my_cmap = copy.copy(plt.get_cmap(cmap))
    my_cmap.set_bad('darkgray')

    # Create list f colors for each area from the ABA color style
    if isinstance(data.index, pd.MultiIndex):
        midIds = data.index.get_level_values('mid').to_list()
    else:
        midIds = data.index.to_list()
    colorList = atlas.ids_to_colors(midIds,color_model='rgb_norm')

    # Create the plot
    axObj = sns.clustermap(data,
                   figsize=figsize,
                   row_colors=colorList,
                   cmap=my_cmap,
                   vmin = vmin,
                   vmax = vmax,
                   mask = data.isna(),
                   row_cluster=False,
                   col_cluster=False,
                   yticklabels=False,
                   dendrogram_ratio=0.01,
                   colors_ratio=0.05,
                   z_score=None,
                   cbar_pos=(1, .3, .03, .3),
                   cbar_kws={"ticks":[vmin,vmax]}
                   )

    # Customize X and Y axes
    axObj.ax_heatmap.set_xlabel('Mice',fontdict={'fontsize':28*fontScaling})
    axObj.ax_heatmap.yaxis.set_label_coords(-0.15,0.5)
    axObj.ax_heatmap.set_xticklabels([])
    axObj.ax_heatmap.set_ylabel('Brain Areas',fontdict={'fontsize':28*fontScaling})
    
    # Change ticks in the colorbar
    axObj.ax_cbar.tick_params(labelsize=20*fontScaling)

    # Add a title if specified
    if title:
        axObj.ax_heatmap.set_title(title, fontsize=28*fontScaling)

    return axObj


def coarseOntologyBarplot(data:pd.DataFrame, atlas:object, cmap:str='PuBu',
    figsize:tuple=(4,6), fontScaling:float=1, xlabel=None, title:str=None, 
    areaNames=True):
    """
    Creates a barplot with mean and SEM as well as single animals.
    You can customize parts of the plot with some of the additional arguments or
    you can use the returned graphical object to customize it

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        cmap:str='PuBu',
            Colormap
        figsize:tuple=(5,12)
            Size of the figure
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure
        xlabel:str=None
            String to use as the X axis title
        title:str=None
            String to use as the plot title
        areaNames:bool=True
            Wether or not to plot area names on the left Y axis

        Returns
        -------
        axObj:obj
            Graphical object normally returned by seaborn.
    """

    # Style. white background and no ticks
    sns.set_style('white')

    # Set up colors for the bars and single animals
    colormap = cm.get_cmap(cmap)
    barColor = colormap(0.4)
    animalColor = colormap(0.85)

    # Create the figure and axes objects
    f, ax = plt.subplots(figsize=figsize)

    # Background Bar Plot
    meltedData = data.melt(ignore_index=False, value_name='metricValue')
    bp = sns.barplot(
        data=meltedData,
        ax=ax,
        y=atlas.ids_to_names(meltedData.index.tolist()),
        x='metricValue',
        orient="h",
        alpha=1,
        linewidth=.5,
        edgecolor="black",
        color=barColor,
        errorbar='se'
    )
    # Plot single animals
    sns.stripplot(
        ax=ax,
        y = atlas.ids_to_names(meltedData.index.tolist()),
        x = 'metricValue',
        data=meltedData,
        size=6,
        orient="h",
        jitter= True,
        alpha=.8,
        linewidth=0.6,
        edgecolor="black",
        color=animalColor,
    )

    # Customize Plot
    bp.xaxis.grid(True)
    maxValue = int(np.ceil(meltedData['metricValue'].max()))
    xticksLocations = [x for x in range(maxValue)]
    ax.xaxis.set_major_locator(FixedLocator(xticksLocations))
    ax.xaxis.set_tick_params(labelsize=22*fontScaling)
    ax.yaxis.set_tick_params(labelsize=18*fontScaling)
    if title:
        ax.set_title(title, fontsize=25*fontScaling)
    if xlabel:
        ax.xaxis.set_label_text(xlabel, fontsize=20*fontScaling)
    if not areaNames:
        ax.yaxis.set_ticks([])
    sns.despine(bottom=True, left=False)

    return ax    


def metricsCorrelation(data:pd.DataFrame, atlas:object, x:str='diffFluo', y:str='energy',
    txtLoc:str='br', ax:object=None, fontScaling:float=1, xlabel:str=None, ylabel:str=None,
    title:str=None):
    """
    Creates a scatter plot with each brain area in "data" by plotting one metric
    against the other.
    You can customize parts of the plot with some of the additional arguments or
    you can use the returned graphical object to customize it

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        x:str='diffFluo'
            Name of the variable in data to plot in the X axis
        y:str='energy'
            Name of the variable in data to plot in the Y axis
        txtLoc:str='br'
            Location of the text of the correlation results.
            'br': bottom-right
            'tl': top-left
        ax:object=None
            axes object where to put the plot
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure
        xlabel:str=None
            String to use as the X axis title
        ylabel:str=None
            String to use as the Y axis title
        title:str=None
            String to use as the plot title

        Returns
        -------
        axObj:obj
            Graphical object normally returned by seaborn.
    """

    # Style. white background and no ticks
    sns.set_style('white')

    if not ax:
        f, ax = plt.subplots(figsize=(4,4))

    # Dictionary that specifies colors for each brain area
    midColorsDict = {id:atlas.ids_to_colors([id],color_model='rgb_norm')[0] for id in atlas.get_midontology_structures_ids()}

    sns.scatterplot(
        ax = ax,
        data=data,
        x=x,
        y=y,
        hue = data.index.get_level_values('mid').tolist(),
        palette=midColorsDict,
        linewidth=0.3,
        edgecolor='k',
        s=150,
        legend=False,
        alpha=1,
    )

    # Setup axis limits
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
    ax.plot(lims, lims, 'darkgray', alpha=0.75, zorder=0)
    ax.set_aspect('equal', adjustable='box')

    # Sprearman Correlation
    coeff, pval = spearmanr(data[x],data[y])
    textColor = cm.get_cmap("OrRd")(0.85)
    lineColor = cm.get_cmap("OrRd")(0.9)

    # Parse text location
    xPos, yPos, align = parseTextLocation(txtLoc)

    # Print text for correlation results
    plt.text(x = xPos, y = yPos+0.12,
        s=r"$r_{s}$" + f"$ = ${coeff:0.2}",
        fontsize=22*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )
    plt.text(x = xPos, y = yPos,
        s=f"$p = ${pval:0.1e}".replace('e-0','e-'),
        fontsize=22*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )

    if pval<0.05:
        huber = HuberRegressor(epsilon=2)
        huber.fit(data[x].values.reshape(-1, 1), data[y].values.ravel())
        newX = np.linspace(data[x].min(), data[x].max(), num=10)
        newY = huber.predict(newX.reshape(-1, 1))
        ax.plot(newX, newY, color=cm.get_cmap('PuBu')(0.9), alpha=0.8, linewidth=3)

    # Axis customization
    ax.set_aspect('equal', adjustable='datalim')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.xaxis.set_tick_params(labelsize=18*fontScaling)
    ax.yaxis.set_tick_params(labelsize=18*fontScaling)
    if xlabel:
        ax.xaxis.set_label_text(xlabel, fontsize=22*fontScaling)
    else:
        ax.xaxis.set_label_text("")
    if ylabel:
        ax.yaxis.set_label_text(ylabel, fontsize=22*fontScaling)
    else:
        ax.yaxis.set_label_text("")

    if title:
        ax.set_title(title, fontsize=26*fontScaling)

    sns.despine()

    return ax


def metricsWithErrors(data:pd.DataFrame, atlas:object, x:str='diffFluo', y:str='energy',
    err_x:str=None, err_y:str=None, ax:object=None, fontScaling:float=1, 
    xlabel:str=None, ylabel:str=None, title:str=None, annotations:bool=True):
    """
    Creates a scatter plot with each major brain area in by plotting one metric
    against the other and adds both vertical and horizontal errorbars.
    You can customize parts of the plot with some of the additional arguments or
    you can use the returned graphical object to customize it

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        x:str='diffFluo'
            Name of the variable in data to plot in the X axis
        y:str='energy'
            Name of the variable in data to plot in the Y axis
        err_x:str=None
            Name of the variable in data to plot errors in the X axis
        err_y:str=None
            Name of the variable in data to plot errors in the Y axis
        ax:object=None
            axes object where to put the plot
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure
        xlabel:str=None
            String to use as the X axis title
        ylabel:str=None
            String to use as the Y axis title
        title:str=None
            String to use as the plot title
        annotations:bool=True
            Whether to plot annotations of area names for each dot

        Returns
        -------
        axObj:obj
            Graphical object normally returned by seaborn.
    """

    # Style. white background and no ticks
    sns.set_style('white')

    if not ax:
        f, ax = plt.subplots(figsize=(4,4))

    # Create a dictionary for ABA-style colors for each coarse area
    coarseIDs = data.index.get_level_values('coarse').unique()
    colors = atlas.ids_to_colors(coarseIDs, color_model='rgb_norm')
    colorDict = {id:color for id,color in zip(coarseIDs,colors)}

    # Main scatterplot
    sns.scatterplot(
        ax = ax,
        data = data,
        x=x,
        y=y,
        hue = data.index.tolist(),
        palette= colorDict,
        linewidth=0.5,
        edgecolor='k',
        s=150,
        legend=False,
        )

    # Plot X and Y errorbars
    plt.errorbar(
        x=data[x],
        y=data[y],
        xerr=data[err_x],
        yerr=data[err_y],
        fmt='None',
        ecolor='dimgrey',
        elinewidth=1.4,
    )

    # Plot the bisector
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
    ax.plot(lims, lims, 'darkgray', alpha=0.75, zorder=0)

    # Axis Labels
    if xlabel:
        ax.xaxis.set_label_text(xlabel, fontsize=22*fontScaling)
    else:
        ax.xaxis.set_label_text("")
    if ylabel:
        ax.yaxis.set_label_text(ylabel, fontsize=22*fontScaling)
    else:
        ax.yaxis.set_label_text("")
    # Title
    if title:
        ax.set_title(title, fontsize=26*fontScaling)

    # Plot annotations
    if annotations:
        for index, row in data.iterrows():
            label = atlas.ids_to_names([index])[0]
            xText = row[x]
            yText = row[y]
            plt.annotate(label,                 # this is the text
                (xText,yText),                  # these are the coordinates to position the label
                textcoords="offset points",     # how to position the text
                xytext=(0,10),                  # distance from text to points (x,y)
                ha='center',                    # horizontal alignment can be left, right or center
                fontsize=18**fontScaling,
                )                    

    # Axis customization
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.xaxis.set_tick_params(labelsize=18*fontScaling)
    ax.yaxis.set_tick_params(labelsize=18*fontScaling)
    ax.set_aspect('equal',adjustable='box')
    sns.despine()

    return ax


def corticalHeatmap(data:pd.DataFrame, atlas:object, cmap:str='PuBu', vmin:float=None, 
    vmax:float=None, figsize:tuple=(10,2.5), fontScaling:float=1, title:str=None):
    """
    Creates a heatmap of cortical areas divided by layers

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        vmin:float=None
            Minimum value for the heatmap
        vmax:float=None
            Maximum value for the heatmap
        cmap:str='PuBu',
            Colormap name
        figsize:tuple=(10,2.5)
            Size of the figure
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure
        title:str=None
            String to use as the plot title

        Returns
        -------
        axObj:obj
            Graphical object normally returned by seaborn.
    """
    
    # Create list of colors for each area from the ABA color style
    listOfIDs = atlas.acronyms_to_ids(data.columns.get_level_values('region').tolist())
    colorList = atlas.ids_to_colors(listOfIDs, color_model='rgb_norm') 
    
    # Set custom color for NaN values
    my_cmap = copy.copy(plt.get_cmap(cmap))
    my_cmap.set_bad('darkgray')

    axObj = sns.clustermap(
        data,
        figsize = figsize,
        vmin = vmin,
        vmax = vmax,
        row_cluster=False,
        col_cluster=False,
        cmap = my_cmap,
        col_colors= colorList,
        mask = data.isna(),
        xticklabels = data.columns.tolist(),
        # square = True,
        cbar_pos=(.98, .4, 0.01, 0.4)
    )

    # Layer labels and ticks
    axObj.ax_heatmap.tick_params(axis='y', 
        left=True, right=False,
        labelleft = True, labelright=False,
        labelrotation=0,
        direction = 'out',
        labelsize=14*fontScaling,
        )
    axObj.ax_heatmap.yaxis.set_label_position('left')
    axObj.ax_heatmap.set_ylabel('Layer', fontsize=16*fontScaling)

    # Remove x ticks from heatmap
    axObj.ax_heatmap.set_xticks([])
    axObj.ax_heatmap.set_xlabel('')
    # Add back x ticks on the color line on top of the heatmap
    axObj.ax_col_colors.xaxis.set_ticks(np.arange(0.5,len(data.columns.tolist())+0.5,1))
    axObj.ax_col_colors.xaxis.set_ticklabels(data.columns.tolist(), fontsize = 12*fontScaling, rotation = 90)
    axObj.ax_col_colors.xaxis.tick_top()

    # Add a title if specified
    if title:
        axObj.ax_heatmap.set_title(title, fontsize=24*fontScaling)

    return axObj


def primaryAreasBarplot(data:pd.DataFrame, metric:str='energy', 
    cmap:str='PuBu', figsize:tuple=(1.5,5), ax:object=None, ylabel:str=None,
    fontScaling:float=1, title:str=None):
    
    # Style. white background and no ticks
    sns.set_style('white')

    if not ax:
        f, ax = plt.subplots(figsize=figsize)

    # Select only the metric to plot
    toPlot = data.xs(metric, axis=1, level='params')
    # Sort to have "primary" before "associative"
    toPlot = toPlot.sort_index(ascending=False)

    # Colors
    cmap = cm.get_cmap(cmap)
    pointColor = cmap(0.85)
    barColor = cmap(0.4)

    paletteBars = {
        'primary':cmap(0.6),
        'associative':cmap(0.2),
    }

    # Plot bars
    sns.barplot(
        ax=ax,
        data=toPlot.T,
        # color=barColor,
        palette=paletteBars,
        errorbar='se',
    )
    # Plot single animals
    ax.plot(
        toPlot,
        linewidth=1,
        color=pointColor,
        marker='o',
        markersize=7,
        alpha=0.5,
    )
    # Axis customization
    ax.set_xlabel("")
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=22*fontScaling)
    if title:
        ax.set_title(title, fontsize=25*fontScaling)
    # Ticks customization
    ax.set_xticklabels(['Primary','Associative'], rotation=45, ha='right', fontsize=18*fontScaling)
    ax.yaxis.set_tick_params(labelsize=18*fontScaling)
    sns.despine()


def primaryAreasLayersBarplot(data:pd.DataFrame, metric:str='energy', 
    cmap:str='PuBu', figsize:tuple=(4,5), ax:object=None, xlabel:str=None,
    fontScaling:float=1, legendTitle:str=None, title:str=None):

    # Style. white background and no ticks
    sns.set_style('white')
    # Create figure and axes
    if not ax:
        f, ax = plt.subplots(figsize=figsize)

    # Melt data in a vertical table
    meltedData = data.xs(metric, axis=1, level='params')\
        .melt(ignore_index=False).reset_index()

    # Define colors for bars and points
    cmap = cm.get_cmap(cmap)
    paletteBars = {
        'primary':cmap(0.6),
        'associative':cmap(0.2),
    }
    palettePoints = {
        'primary':cmap(0.85),
        'associative':cmap(0.85),
    }

    # Plot bars
    sns.barplot(
        data=meltedData,
        ax=ax,
        orient='h',
        y='layer',
        x='value',
        hue='sensory',
        hue_order=['primary','associative'],
        palette=paletteBars,
        # saturation=.7,
        errorbar='se',
    )
    # Plot single animals
    sns.stripplot(
        data=meltedData,
        ax=ax,
        y='layer',
        x='value',
        hue='sensory',
        hue_order=['primary','associative'],
        palette=palettePoints,
        size=5,
        dodge=True,
        alpha=0.75,
    )

    # Customize legend
    h,l = ax.get_legend_handles_labels()
    if not legendTitle:
        legendTitle='Sensory areas'
    ax.legend(h[2::], l[2::], title=legendTitle, fontsize=14*fontScaling, 
        title_fontsize=16*fontScaling,loc='best', frameon=False)

    # Customize axis
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=20*fontScaling)
    if title:
        ax.set_title(title, fontsize=25*fontScaling)
    ax.set_ylabel("Cortical Layer",fontsize=22*fontScaling)
    ax.tick_params(labelsize=18*fontScaling)

    sns.despine()


def volcano_plot(
                    stat_dataframe:pd.DataFrame,
                    atlas:object,
                    title:str = None,
                    fontScaling:float = 1,
                    figsize:tuple =(7,6)):

    sns.set_style('white')
    colorsDict = {id:atlas.ids_to_colors([id],color_model='rgb_norm')[0] for id in stat_dataframe.index.tolist()}
    _, ax = plt.subplots(figsize = figsize)
    g = sns.scatterplot(   
                    data=stat_dataframe, 
                    x = 'log2FC',
                    y = '-log10p',
                    hue = stat_dataframe.index.tolist(),
                    palette= colorsDict,
                    legend = False, 
                    s = 95, 
                    linewidth=0.25,
                    edgecolor='k', 
                    alpha = 1
                    )
    
    
    if title:
        ax.set_title(title, fontsize=25*fontScaling)
    ax.set_xlabel(r"$\log_{2}{FC}}$", fontsize=20*fontScaling)
    ax.set_ylabel(r"$-\log_{10}{p}}$",fontsize=22*fontScaling)
    ax.tick_params(labelsize=16*fontScaling)


    sns.despine()
    
    return g


def ora_plot(stat_df:pd.DataFrame, 
            hue = 'p_value',
            size = 'percentage_of_term',
            custom_thr = 0.05,
            figsize:tuple = (7, 5),
            title = None,
            ylabel:str = None,
            dotsizes:tuple = (1, 300),
            fontScaling = 1,
            cmap = 'coolwarm_r'):

    sns.set_style('white')


    stat_df = stat_df.copy().reset_index()


    #colormap and main figure
    cmap = cm.get_cmap(cmap)
    fig, axs = plt.subplots(ncols=2, nrows=2, gridspec_kw={'height_ratios':[3,2],'width_ratios': [15,1]}, figsize = figsize)
    gs = axs[0,0].get_gridspec()
    for ax in axs[0:, 0]:
        ax.remove()
    fig.add_subplot(gs[0:, 0])
    axes = fig.get_axes()

    norm = colors.TwoSlopeNorm(vcenter = custom_thr, vmin = 0, vmax=1)


    pal = get_custom_palette(norm,cmap ,stat_df[hue])
    sns.scatterplot(
                    data = stat_df,
                    x = 'enrichment_ratio',
                    y = 'set_name',
                    hue = hue,
                    size= size,
                    edgecolor = 'k',
                    linewidth = .4,
                    palette= pal,
                    ax = axes[2],
                    sizes=dotsizes,
                    zorder =20
                    )
    xlims = axes[2].get_xlim()
    for ind, area in stat_df.iterrows():
        sns.lineplot(y=[ind, ind], x= [xlims[0],area['enrichment_ratio']],
            zorder = -20,
            color = pal[area[hue]],
            linewidth = 3.5,
            ax = axes[2])
    
    axes[2].set_xlim(xlims)
    axes[2].set_title(title, fontsize = 18*fontScaling, y = 1.09)
    axes[2].tick_params(axis = 'both',labelsize = 14*fontScaling)
    axes[2].set_xlabel('Enrichment ratio', fontsize = 15*fontScaling)
    if ylabel:
        axes[2].set_ylabel(ylabel, fontsize = 15*fontScaling)
    else:
        axes[2].set_ylabel('', fontsize = 15*fontScaling)


    #legend
    h, l = axes[2].get_legend_handles_labels()
    axes[2].legend([],[], frameon=False)
    idx = l.index(size)
    handles = h[idx:]
    labels = []
    labels.append(r'Percentage of set')
    labels.extend([f'{lb:.3}'+r' %' for lb in l[idx+1:]])
    
    axes[0].legend(handles,labels,bbox_to_anchor=(0,1),  loc='upper left',borderaxespad=0)

    axes[0].axis('off')

    # colorbar
    cbar = plt.colorbar(
                mappable=  cm.ScalarMappable( cmap=cmap),
                ticks = [0, .5, 1],
                cax=axes[1])
    if hue == 'p_value':
        title = 'p-value'
    elif hue == 'FDR':
        title = 'FDR'
    axes[1].set_title(title)
    cbar.set_ticklabels([0, custom_thr, 1])

    # fig.tight_layout()

    sns.despine()
    return fig
    

def pca_plot(pca_df:pd.DataFrame, var:tp.Sequence,
                # style = 'sex',

                size = None,
                weightrange:tuple = (22,44),
                dotsizes = (10,400),
                variable = 'sex',
                cmap = 'magma',
                fontScaling = 1,
                title:str = None
                ):
    sns.set_style('white')


    if len(variable) == 1:

        palette = sns.color_palette(cmap)
        colorsCh = [palette[1],palette[5]]



        ax = sns.scatterplot(data=pca_df,
                        x = 'PC1',
                        y= 'PC2',
                        hue = variable[0],
                        size = size,
                        # style = variable[0],
                        size_norm=weightrange,
                        sizes = dotsizes,
                        palette=colorsCh,
                        s = 200
                        )


        for i, areacategory in enumerate(pca_df[variable[0]].unique()):
            set = pca_df.loc[pca_df[variable[0]]==areacategory, ['PC1', 'PC2']]
            points = set.values
            hull = ConvexHull(points)

            p = Polygon(points[hull.vertices,:], alpha=0.3, zorder=-20, 
            facecolor=colorsCh[i]
            )
            ax.add_patch(p)
    elif len(variable)==2:
        s_treat = [s+'-'+t for s, t in zip(pca_df[variable[0]], pca_df[variable[1]])]
        colstr = variable[0]+'-'+variable[1]
        pca_df[colstr] = s_treat

        palette = sns.color_palette(cmap)
        colorsCh = [palette[0],palette[2],palette[4],palette[5]]

        ax = sns.scatterplot(data=pca_df,
                x = 'PC1',
                y= 'PC2',
                hue = colstr,
                size = size,
                style = variable[0],
                size_norm=weightrange,
                sizes = dotsizes,
                palette=colorsCh,
                s = 200
                )
        for i, areacategory in enumerate(pca_df[colstr].unique()):
            set = pca_df.loc[pca_df[colstr]==areacategory, ['PC1', 'PC2']]
            points = set.values
            hull = ConvexHull(points)

            p = Polygon(points[hull.vertices,:], alpha=0.3, zorder=-20, 
            facecolor=colorsCh[i])
            ax.add_patch(p)


    ax = plt.gca()

    ax.set_aspect('equal', adjustable = 'datalim')
    ax.set_xlabel(f"PC1 ({var[0]*100:2.0f}%)", fontdict={'size': 15*fontScaling})
    ax.set_ylabel(f"PC2 ({var[1]*100:2.0f}%)", fontdict={'size': 15*fontScaling})
    if title:
        ax.set_title(title, fontdict={'size': 20*fontScaling}, y = 1.05)
    ax.set_aspect('equal', adjustable = 'box')
    # plt.axis('equal')
    ax.legend(bbox_to_anchor=(1.05,1),  loc='upper left',borderaxespad=0)
    sns.despine()

    return ax


def two_groups_heatmap(
                data:pd.DataFrame,
                atlas:object,
                cmap:str='OrRd',
                col_cluster:bool=True,
                metric:str = 'correlation',
                method:str='average',
                center:float=1,
                vmax:float = 4,
                vmin:float = 0,
                figsize:tuple = (7, 10), 
                fontScaling:float = 1,
                title:str = None
            ):
    '''Heatmap for the visualization of differences between two experimental groups.

    Parameters
    ----------
    data:pd.DataFrame,
    atlas:object,
    cmap:str,
    col_cluster:bool,
    metric:str,
    method:str,
    center:float,
    vmax:float,
    vmin:float = -2,
    figsize:tuple = (7, 10), 
    fontScaling:float = 1,
    title:str = None



    Returns
    -------
    axObj:obj
        Graphical object normally returned by seaborn.
    '''
    sns.set_style('white')

    # Colors for the diferent mice
    treatments = [x[0] for x in data.columns.to_list()]
    columnList = [sns.color_palette('deep')[0] if x=='CTR' else sns.color_palette('deep')[1] for x in treatments]
    colorList = atlas.ids_to_colors(data.index.get_level_values('mid').to_list(), color_model='rgb_norm')

    if col_cluster:
        data = data.dropna()
        mask = None
    else:
        mask = data.isna()
    g = sns.clustermap(
                    data,
                    figsize=figsize,
                    row_cluster=False,
                    col_cluster=col_cluster,
                    cmap=cmap,
                    # center=center,
                    mask = mask,
                    metric=metric,
                    method=method,
                    # z_score=1,
                    dendrogram_ratio=(.12),
                    cbar_pos=(1, .3, .03, .3),
                    row_colors=colorList,
                    col_colors=columnList,
                    vmax = vmax,
                    vmin = vmin, 
                    yticklabels=False
                    )

    ax = g.ax_heatmap
    if col_cluster ==  False:
        ax1 = g.ax_col_colors
        ax1.set_title(title, fontdict={'fontsize':28 * fontScaling})
    else:
        ax1 = g.ax_col_dendrogram
        ax1.set_title(title, fontdict={'fontsize':28* fontScaling})
    
    ax.set_xlabel('Treatment - Mouse',fontdict={'fontsize':25* fontScaling})

    ax.set_ylabel('Brain Areas',fontdict={'fontsize':25* fontScaling})
    ax.tick_params(labelsize = 18*fontScaling)
    g.ax_cbar.tick_params(labelsize = 15*fontScaling)
    ax.yaxis.set_label_coords(-.13,0.5)
    ax.xaxis.set_label_coords(0.5, -0.18)
    return g


def get_custom_palette(norm:colors.Normalize, cmap:str, vals):
    
    cmap_obj = plt.get_cmap(cmap)
    pal = {v:cmap_obj(norm(v)) for v in vals}
    return pal 


def connectomeCorrelationScatterplot(data:pd.DataFrame, atlas:object, layer:str='4', txtLoc:str='br', 
    ax:object=None, fontScaling:float=1, xlabel:str=None, ylabel:str=None, title:str=None):
    """
    Creates a scatter plot with each brain area in "data" by plotting the data
    specified by the column name "layer" versus the column "afferent".
    You can customize parts of the plot with some of the additional arguments or
    you can use the returned graphical object to customize it

        Parameters
        ----------
        data:pd.DataFrame
            Dataframe of a single metric of interest
        atlas:object
            Instance of the Atlas class defined in AbaTool.py
        layer:str='4'
            Name of the variable in data to plot in the X axis ('1','2/3','4','5','6')
        txtLoc:str='br'
            Location of the text of the correlation results.
            'br': bottom-right
            'tl': top-left
        ax:object=None
            axes object where to put the plot
        fontScaling:float=1
            Multiplicative factor for all the fonts in the figure
        xlabel:str=None
            String to use as the X axis title
        ylabel:str=None
            String to use as the Y axis title
        title:str=None
            String to use as the plot title

        Returns
        -------
        None
    """

    # Style. white background and no ticks
    sns.set_style('white')

    if not ax:
        f, ax = plt.subplots(figsize=(4,4))

    # Dictionary that specifies colors for each brain area
    areaIds = data.index.get_level_values('mid')
    colorsDict = {id:atlas.ids_to_colors([id],color_model='rgb_norm')[0] for id in areaIds}

    sns.scatterplot(
        data=data,
        ax=ax,
        x='afferents',
        y=layer,
        s=150,
        hue=data.index.get_level_values('mid'),
        palette=colorsDict,
        edgecolor='black',
        linewidth=0.4,
        legend=False,
    )

    # Sprearman Correlation
    # coeff, pval = spearmanr(data['afferents'], data[layer], nan_policy='omit')
    coeff, pval = pearsonr(data['afferents'], data[layer])
    textColor = cm.get_cmap("OrRd")(0.85)
    lineColor = cm.get_cmap("OrRd")(0.9)
    # Parse text location
    if txtLoc=='br':
        xPos = 1
        yPos = 0.05
        align = 'right'
    elif txtLoc=='tl':
        xPos = 0.05
        yPos = 0.80
        align = 'left'
    # Print text for correlation results
    plt.text(x = xPos, y = yPos+0.12,
        s=r"$r$" + f"$ = ${coeff:0.2}",
        fontsize=20*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )
    plt.text(x = xPos, y = yPos,
        s=f"$p = ${pval:0.1e}".replace('e-0','e-'),
        fontsize=20*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )

    if pval<0.05:
        huber = HuberRegressor(epsilon=2)
        huber.fit(data['afferents'].values.reshape(-1, 1), data[layer].values.ravel())
        newX = np.linspace(data['afferents'].min(), data['afferents'].max(), num=10)
        newY = huber.predict(newX.reshape(-1, 1))
        ax.plot(newX, newY, color=cm.get_cmap('PuBu')(0.9), alpha=0.8, linewidth=3)

    # Axis customization
    ax.set_aspect('equal', adjustable='datalim')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.xaxis.set_tick_params(labelsize=16*fontScaling)
    ax.yaxis.set_tick_params(labelsize=16*fontScaling)
    if xlabel:
        ax.xaxis.set_label_text(xlabel, fontsize=20*fontScaling)
    else:
        ax.xaxis.set_label_text("")
    if ylabel:
        ax.yaxis.set_label_text(ylabel, fontsize=20*fontScaling)
    else:
        ax.yaxis.set_label_text("")

    if title:
        ax.set_title(title, fontsize=26*fontScaling)

    sns.despine()


def correlationWithGene(atlas:object, x:pd.Series=None, y:pd.Series=None, txtLoc:str='tl', 
    pval:float=None, corr_spearman:float=None,ax:object=None, fitLine:bool=True, 
    fontScaling:float=1, xlabel:str=None, ylabel:str=None, title:str=None):

    sns.set_style('white')

    if not ax:
        f,ax = plt.subplots(figsize=(4,4))

    # Remove elements == 0 since they give inf when log10
    df = pd.concat([x,y], axis=1)
    df.columns = ['x','y']
    df = df.loc[(df!=0).all(axis=1)]
    x = df['x']; y=df['y']
    
    # Calculate the Log10 of the variables
    xlog = np.log10(x).values
    ylog = np.log10(y).values

    ax.set_xscale('log')
    ax.set_yscale('log')

    # Dictionary that specifies colors for each brain area
    midColorsDict = {id:atlas.ids_to_colors([id],color_model='rgb_norm')[0] for id in atlas.get_midontology_structures_ids()}

    # Plot points
    sns.scatterplot(
        ax=ax,
        x=x,
        y=y,
        s=25,
        palette=midColorsDict,
        hue=x.index,
        edgecolor='dimgrey',
        legend=False,
    )

    # Adjust axes and get the limits
    ax.set_aspect('equal', adjustable='datalim')
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()

    # KDE plot
    sns.kdeplot(
        x=x,
        y=y,
        fill=False,
        cmap='PuBu',
        levels=4,
        thresh=0.1,
        log_scale=True,
        zorder=-2,
        norm=colors.Normalize(vmin=-.5, vmax=.6),
        ax=ax
    )  

    # Sprearman Correlation
    xPos, yPos, align = parseTextLocation(txtLoc)
    textColor = cm.get_cmap("OrRd")(0.85)
    # Print text for correlation results
    plt.text(x = xPos, y = yPos+0.12,
        s=r"$r_{s}$" + f"$ = ${corr_spearman:0.2}",
        fontsize=22*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )
    plt.text(x = xPos, y = yPos,
        s=f"$p = ${pval:0.1e}".replace('e-0','e-'),
        fontsize=22*fontScaling,
        transform=ax.transAxes,
        horizontalalignment=align,
        color = textColor if pval <0.05 else 'dimgray'
    )

    # If to plot a regression line or not
    if fitLine:
        huber = HuberRegressor(epsilon=1.35)
        huber.fit(xlog.reshape(-1, 1), ylog.ravel())
        newX = np.linspace(np.min(xlog), np.max(xlog), num=3)
        newY = huber.predict(newX.reshape(-1, 1))
        plt.plot(10**newX,10**newY, color=cm.PuBu(0.9),alpha=0.8, linewidth=3)

    # Set axis limits back after KDE plot
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    # Axis customization
    ax.xaxis.set_tick_params(labelsize=14*fontScaling)
    ax.yaxis.set_tick_params(labelsize=14*fontScaling)
    if xlabel:
        ax.xaxis.set_label_text(xlabel, fontsize=22*fontScaling)
    else:
        ax.xaxis.set_label_text("")
    if ylabel:
        ax.yaxis.set_label_text(ylabel, fontsize=22*fontScaling)
    else:
        ax.yaxis.set_label_text("")

    if title:
        ax.set_title(title, fontsize=26*fontScaling, style='italic')

    sns.despine()


def colocalizationDoubleBarplot(data:pd.DataFrame, atlas:object, x_left:str='pvPositive_pnn',
    x_right:str='wfaPositive_pv', figsize:tuple=(6,6), cmaps:list=['bone_r', 'pink_r'],
    xlabel_left:str=None, xlabel_right:str=None, title_left:str=None,
    title_right:str=None, fontScaling:float=1, adaptiveHeight:bool=False, dots:bool=True):
    
    # Style. white background and no ticks
    sns.set_style('white')
    # Create figure and axes
    if adaptiveHeight:
        figsize = (figsize[0], data.shape[0]*0.35)
    f, axs = plt.subplots(ncols= 2, sharey=True, figsize=figsize)

    cmaps = [cm.get_cmap(cmaps[0]),cm.get_cmap(cmaps[1])]

    # Plot the left and right barplots
    for i, ax in enumerate(axs):

        # Select which variable to plot
        varToPlot = x_left if i==0 else x_right

        # Select color for this side of the graph
        barColor = cmaps[i](0.5)
        animalColor = cmaps[i](0.9)

        # Plot background gray bar
        sns.barplot(
            ax=ax,
            y = atlas.ids_to_names(data.index.tolist()),
            x = [100] * data.shape[0],       # All bars go up to 100
            orient="h",
            alpha=1,
            linewidth=0,
            facecolor='gainsboro',
        )

        meltedData = data.xs(varToPlot, axis='columns', level='params').melt(ignore_index=False)
        # Plot experimental data bars
        sns.barplot(
            ax=ax,
            y = atlas.ids_to_names(meltedData.index.tolist()),
            x = 'value',
            data = meltedData,
            orient = "h",
            saturation = 1,
            alpha=1,
            linewidth=0,
            color=barColor,
            errorbar='se'
        )
        # Single animal dots
        if dots:
            sns.stripplot(
                ax=ax,
                y = atlas.ids_to_names(meltedData.index.tolist()),
                x = 'value',
                data=meltedData,
                size=5,
                orient="h",
                jitter= True,
                alpha=1,
                linewidth=0.3,
                edgecolor="black",
                color=animalColor,
            )
        ax.bar_label(
            ax.containers[1],
            fmt='%.1f',
            label_type='edge',
            padding=1,
            fontsize=13*fontScaling,
        )

    # Flip the X orientation opf the left plot
    axs[0].invert_xaxis() 
    # Make the plot adjacent
    plt.subplots_adjust(wspace=0)

    sns.despine(bottom=False, left=True,  ax = axs[0])
    sns.despine(bottom=False, left=False, ax = axs[1])

    # Customize ticks
    axs[0].xaxis.set_tick_params(labelsize=16*fontScaling)
    axs[0].yaxis.set_tick_params(labelsize=18*fontScaling)
    axs[1].xaxis.set_tick_params(labelsize=16*fontScaling)

    # Customize axes
    if xlabel_left:
        axs[0].xaxis.set_label_text(xlabel_left, fontsize=18*fontScaling)
    else:
        axs[0].xaxis.set_label_text("")
    if xlabel_right:
        axs[1].xaxis.set_label_text(xlabel_right, fontsize=18*fontScaling)
    else:
        axs[1].xaxis.set_label_text("")

    if title_left:
        axs[0].set_title(title_left, fontsize=20*fontScaling)
    if title_right:
        axs[1].set_title(title_right, fontsize=20*fontScaling)


def energyDiffuseDoubleBarplot(data:pd.DataFrame, atlas:object, x_left:str='pvPositive_pnn',
    x_right:str='wfaPositive_pv', figsize:tuple=(6,6), cmaps:list=['bone_r', 'pink_r'],
    xlabel_left:str=None, xlabel_right:str=None, title:str = None, fontScaling:float=1, adaptiveHeight:bool=False, dots:bool=True):
    
    # Style. white background and no ticks
    sns.set_style('white')
    # Create figure and axes
    if adaptiveHeight:
        figsize = (figsize[0], data.shape[0]*0.35)
    f, axs = plt.subplots(ncols= 2, sharey=True, figsize=figsize)

    cmaps = [cm.get_cmap(cmaps[0]),cm.get_cmap(cmaps[1])]

    # Plot the left and right barplots
    for i, ax in enumerate(axs):

        # Select which variable to plot
        varToPlot = x_left if i==0 else x_right

        # Select color for this side of the graph
        barColor = cmaps[i](0.5)
        animalColor = cmaps[i](0.9)

        meltedData = data.xs(varToPlot, axis='columns', level='params').melt(ignore_index=False)


        bp  = sns.barplot(
                data=meltedData.dropna(),
                ax=ax,
                y=atlas.ids_to_names(meltedData.dropna().index.tolist()),
                x='value',
                orient="h",
                alpha=1,
                linewidth=.5,
                edgecolor="black",
                color=barColor,
                errorbar='se'
            )
        if dots:
            # Plot single animals
            sns.stripplot(
                    ax=ax,
                    y = atlas.ids_to_names(meltedData.dropna().index.tolist()),
                    x = 'value',
                    data=meltedData.dropna(),
                    size=6,
                    orient="h",
                    jitter= True,
                    alpha=.8,
                    linewidth=0.6,
                    edgecolor="black",
                    color=animalColor,
                )
        bp.xaxis.grid(True)



    # Flip the X orientation opf the left plot
    axs[0].invert_xaxis() 
    # Make the plot adjacent
    plt.subplots_adjust(wspace=0)
    
    sns.despine(bottom=False, left=True,  ax = axs[0])
    sns.despine(bottom=False, left=False, ax = axs[1])

    # Customize ticks
    axs[0].xaxis.set_tick_params(labelsize=16*fontScaling)
    axs[0].yaxis.set_tick_params(labelsize=18*fontScaling)
    axs[1].xaxis.set_tick_params(labelsize=16*fontScaling)

    # Customize axes
    if xlabel_left:
        axs[0].xaxis.set_label_text(xlabel_left, fontsize=18*fontScaling)
    else:
        axs[0].xaxis.set_label_text("")
    if xlabel_right:
        axs[1].xaxis.set_label_text(xlabel_right, fontsize=18*fontScaling)
    else:
        axs[1].xaxis.set_label_text("")


def colocProbabilityPlot(data:pd.DataFrame, ax:object=None, xlabel:str=None, 
    ylabel:str=None, txtLoc:str='tl', fontScaling:float=1, cmap:str='PuBu',
    miceColor=None, avgColor=None, ylim:tuple=None, xticks:bool=True, yticks:bool=True, 
    title:str=None, singleMice:bool=True, txtStat:bool=True):

    data = data.copy()
    if not miceColor:
        cmap = cm.get_cmap(cmap)
        miceColor = cmap(0.3)
    if not avgColor:
        cmap = cm.get_cmap(cmap)
        avgColor = cmap(0.9)

    if not ax:
        f, ax = plt.subplots(figsize=(3,4))

    if singleMice:
        ax.plot(
            data.T,
            linewidth=1,
            color = miceColor,
            zorder=-1
        )

    melted = data.melt(ignore_index=False)
    sns.pointplot(
        data=melted,
        ax=ax,
        x='intClass',
        y='value',
        errorbar='se',
        color= avgColor,
        scale=1.1,
    )

    # Do the linear fit
    x = list(range(1,data.shape[1]+1,1))
    y = data.mean()
    deg = 1
    fitCoefficients = np.polyfit(x, y, deg)

    # Do the ANOVA
    toAnalize = data.dropna(axis=0, how='any')
    model = anova.AnovaRM(
        data=toAnalize.melt(ignore_index=False).reset_index(),
        depvar='value',
        subject='mouse',
        within=['intClass']
    ).fit()
    t = model.anova_table.to_numpy()
    DF1 = t[0,1]
    DF2 = t[0,2]
    F = t[0,0]
    pval = t[0,3]

    xPos, yPos, align = parseTextLocation(txtLoc)

    if txtStat:
        signifColor = cm.get_cmap("OrRd")(0.9)
        plt.text(x = xPos, y = yPos,
            s=f"F$_{{({DF1:.0f},{DF2:.0f})}}$={F:.2f}",
            fontsize=17*fontScaling,
            transform=ax.transAxes,
            horizontalalignment=align,
            color = signifColor if pval <0.05 else 'dimgray'
        )
        plt.text(x = xPos, y = yPos+0.1,
            s=f"p ={pval:.3f}" if pval>=0.001 else "P<0.001",
            fontsize=17*fontScaling,
            transform=ax.transAxes,
            horizontalalignment=align,
            color = signifColor if pval <0.05 else 'dimgray'
        )
        plt.text(x = xPos, y = 0.05,
            s=f"Y = {fitCoefficients[0]:.2f}*X + {fitCoefficients[1]:.2f}",
            fontsize=16*fontScaling,
            transform=ax.transAxes,
            horizontalalignment=align,
            color = 'dimgray'
        )


    # Customize axis
    sns.despine()
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=18*fontScaling)
    else:
        ax.set_xlabel("", fontsize=18*fontScaling)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=18*fontScaling)
    else:
        ax.set_ylabel("", fontsize=18*fontScaling)
    if ylim:
        ax.set_ylim(ylim)
    if not xticks:
        ax.set_xticks([])
    if not yticks:
        ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=20*fontScaling)
        
    ax.tick_params(labelsize=16*fontScaling)


def parseTextLocation(txtLocation:str='tl'):
    
    # Parse text location
    if txtLocation=='br':
        xPos = 1
        yPos = 0.05
        align = 'right'
    elif txtLocation=='tl':
        xPos = 0.05
        yPos = 0.80
        align = 'left'

    return xPos, yPos, align

# ------------------------------------------------------------------------------
# ACCESSORY FUNCTIONS
# ------------------------------------------------------------------------------
COLORMAP_FILE = ''

class ColorMap2D:
    def __init__(self, filename: str, transpose=False, reverse_x=False, reverse_y=False, xclip=None, yclip=None):
        """
        Maps two 2D array to an RGB color space based on a given reference image.
        Args:
            filename (str): reference image to read the x-y colors from
            rotate (bool): if True, transpose the reference image (swap x and y axes)
            reverse_x (bool): if True, reverse the x scale on the reference
            reverse_y (bool): if True, reverse the y scale on the reference
            xclip (tuple): clip the image to this portion on the x scale; (0,1) is the whole image
            yclip  (tuple): clip the image to this portion on the y scale; (0,1) is the whole image
        """
        self._colormap_file = filename or COLORMAP_FILE
        self._img = plt.imread(self._colormap_file)
        self._img = self._img[:,:,0:3]
        if transpose:
            self._img = self._img.transpose()
        if reverse_x:
            self._img = self._img[::-1,:,:]
        if reverse_y:
            self._img = self._img[:,::-1,:]
        if xclip is not None:
            imin, imax = map(lambda x: int(self._img.shape[0] * x), xclip)
            self._img = self._img[imin:imax,:,:]
        if yclip is not None:
            imin, imax = map(lambda x: int(self._img.shape[1] * x), yclip)
            self._img = self._img[:,imin:imax,:]
        if issubclass(self._img.dtype.type, np.integer):
            self._img = self._img / 255.0

        self._width = len(self._img)
        self._height = len(self._img[0])

        self._range_x = (0, 1)
        self._range_y = (0, 1)


    @staticmethod
    def _scale_to_range(u: np.ndarray, u_min: float, u_max: float) -> np.ndarray:
        return (u - u_min) / (u_max - u_min)

    def _map_to_x(self, val: np.ndarray) -> np.ndarray:
        xmin, xmax = self._range_x
        val = self._scale_to_range(val, xmin, xmax)
        rescaled = (val * (self._width - 1))
        return rescaled.astype(int)

    def _map_to_y(self, val: np.ndarray) -> np.ndarray:
        ymin, ymax = self._range_y
        val = self._scale_to_range(val, ymin, ymax)
        rescaled = (val * (self._height - 1))
        return rescaled.astype(int)

    def __call__(self, val_x, val_y):
        """
        Take val_x and val_y, and associate the RGB values 
        from the reference picture to each item. val_x and val_y 
        must have the same shape.
        """
        if val_x.shape != val_y.shape:
            raise ValueError(f'x and y array must have the same shape, but have {val_x.shape} and {val_y.shape}.')
        self._range_x = (np.amin(val_x), np.amax(val_x))
        self._range_y = (np.amin(val_y), np.amax(val_y))
        x_indices = self._map_to_x(val_x)
        y_indices = self._map_to_y(val_y)
        self.i_xy = np.stack((x_indices, y_indices), axis=-1)
        self.rgb = np.zeros((*val_x.shape, 3))
        for indices in np.ndindex(val_x.shape):
            self.img_indices = tuple(self.i_xy[indices])
            self.rgb[indices] = self._img[self.img_indices]
        return self.rgb

    def generate_cbar(self, nx=100, ny=100):
        "generate an image that can be used as a 2D colorbar"
        x = np.linspace(0, 1, nx)
        y = np.linspace(0, 1, ny)
        return self.__call__(*np.meshgrid(x, y))