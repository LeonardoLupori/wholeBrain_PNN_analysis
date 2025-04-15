[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7991607.svg)](https://doi.org/10.5281/zenodo.7991607)
# wholeBrain_PNN_analysis

Analysis code to reproduce all figures of the paper [_A Comprehensive Atlas of Perineuronal Net Distribution and Colocalization with Parvalbumin in the Adult Mouse Brain_](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00799-4)

## Dependencies

To run the scripts in this repository you will need:

1. To download preprocessed data
2. To install the correct version of Python and additional packages

### Download preprocessed data

In order to run the analysis code you will need to download a set of preprocessed data. You can find the data files in the latest release under the _assets_ menu at this [link](https://github.com/LeonardoLupori/wholeBrain_PNN_analysis/releases/latest)

> **Please Note**  
> If you want to download all the **raw** microscopy data that are related to this paper, you can download the entire dataset on Zenodo [HERE](https://doi.org/10.5281/zenodo.7419283).

### Python Environment

To run the scripts in this repository you will need python 3.8 and the correct version of several additional packages. You can easily install everything in a new virtual environment using anaconda.

To do this, use the following commands:

```bash
conda create -n pnnWholeBrain python==3.8
```

Then activate the environment and install all the required packages

```bash
conda activate pnnWholeBrain
pip install -r requirements.txt
```

## The paths.ini file

The repository contains a `paths.ini` configuration file, which is used to parse user-specified paths to the various functions in the notebooks, making collaborative programming easier. These paths are specified in the `paths.ini` file and read by the pathParser object in the dataIO module. They are then assigned to dedicated variables. Alternatively, paths can be assigned directly to these variables rather than being specified in the `paths.ini` file.

To run the notebooks download the `DATA` folder from the assets of this repo and manually put the full path of all the requested files and folders in the `paths.ini` file.

## wholebrain_tools package

Analysis in this work was conducted with the help of custom code wrapped in a python package. This package contains the following modules:

- **aba**: handling of the Allen Brain Atlas
- **stats**: tools for statistical analysis
- **graphics**: tools for custom plots
- **dataIO**: file handling
- **background**: functions to analyze the control dataset
- **pvnegative**: function to analyze the energy of a subset of PNNs (PV<sup>-</sup>)

## Figure 2

### Distribution of WFA-positive PNNs throughout the entire mouse brain

Plots in this figures provide a graphical visualization of the distriburion of WFA-positive perineuronal nets throughout the mouse brain.

The folder contains the following notebooks:

- **figure_02_mainVisualizations.ipynb** all the plots of the figure (scatterplots, heatmaps), representing the relationship of WFA diffuse staining and PNN energy as well as raw data of individual animals.
- **figure_02_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_02_brainRenders.ipynb. Saves a .csv file used for plotting by the brainRenders notebook
- **figure_02_brainRenders.ipynb** heatmaps representing on brain coronal sections WFA diffuse fluorescence and PNN energy at mid-ontology resolution.

## Figure 3

### Brain-wide interactions between PNNs and PV cells

Plots in this figures provide a graphical visualization of the colocalization of PNNs and PV cells throughout the mouse brain.

The folder contains the following notebooks:

- **figure_03_colocalization.ipynb** barplots to visualize PNN/PV colocalization at the coarse level
- **figure_03_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_03_brainRenders.ipynb. Saves a .csv file used for plotting by the brainRenders notebook
- **figure_03_colocalizationBrainRenders.ipynb** heatmaps representing on brain coronal sections the percentage of PNNs surrounding a PV cells and the percentage of PV cells ensheated by PNNs.
- **figure_03_correlation_pv_wfa.ipynb** scatterplots representing the relationship of WFA diffuse staining/PNN energy and PV energy.

## Figure 4

### PNNs aggregation depends on PV expression levels

Plots in this figures explore the relationship between probability of having a PNN and PV expression levels throughout the mouse brain.

The folder contains the following notebooks:

- **figure_04_colocProbability.ipynb**  scatterplots representing the relationship between PNN expression and intensity staining of individual PV cells.

## Figure 5

### Principles of PNN organization in cortical areas

Plots in this figures explore the expression of PNNs in the cortex.

The folder contains the following notebooks:

- **figure_05_corticalHeatmap.ipynb** heatmaps representing PNN expression in cortical areas.
- **figure_05_primarySecondary.ipynb** barplots representing PNN expression in primary vs associative sensory areas.
- **figure_05_functionalRegions.ipynb** scatterplots representing the relationship of PNN expression and intensity staining of individual PV cells.
- **figure_05_connectome.ipynb** scatterplots representing the relationship between PNN expression in layers of the sensory cortex and density of thalamic afference.

## Figure 6

### Gene expression correlates of PNNs

Plots in this figures explore correlation of PNN energy with known PNN molecular markers and the biological processes in which the genes correlated or anticorrelated with WFA staining are involved.

- **figure_06_downloadGeneExprABA.ipynb** download of the AGEA dataset as a .csv file (see [Lein et al., 2007](https://doi.org/10.1038/nature05453)).
- **figure_06_geneOntology.ipynb** production of gene lists for gene ontology analysis and plotting of the results.
- **figure_06_matrisome.ipynb** production of gene lists for overrepresentation analysis in the matrisome dataset [Naba et al., 2016](https://doi.org/10.1016/j.matbio.2015.06.003).
- **figure_06_plotMarkers.ipynb** scatterplots representing PNN/PV levels as a function of the expresssion of known molecular markers.
- **figure_06_preprocess_dataWFA.ipynb** correlation analysis of WFA fluorescence/PNN energy and gene expression data of AGEA dataset.
- **figure_06_preprocess_dataPV.ipynb** correlation analysis of PV energy and AGEA dataset.

## Figure S01

### Evaluation of the scoring model effectiveness and validation of the defined staining metrics

Folder content:

#### The scores assigned by the scoring models correlate with raters’ agreement

Plots showing the performances of our scoring model.

- **figure_S01_networkMetrics.ipynb** boxplots representing the score assigned by the model to cells at different levels of agreement.

####  Negligible contribution of WFA background fluorescence to the WFA Diffuse Fluorescence

Plots showing a comparison between WFA-stained slices and negative control slices that lacked the use of biotinylated WFA (but not the streptavidin-conjugated fluorescent probe) in the staining protocol.

- **figure_S01_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_S02_brainRenders.ipynb. Saves a .csv file used for plotting by the brainRenders notebook.
- **figure_S01_brainRenders.ipynb** heatmaps representing on brain coronal sections average pixel intensity at mid-ontology resolution.

## Figure S02

### PNN energy and WFA diffuse fluorescence measurements for medium-resolution brain areas grouped by their major subdivision

Folder content:

- **figure_S02.ipynb** Barplots showing PNN energy and WFA diffuse fluorescence at mid-ontology level.

## Figure S03

### Distribution of PV-positive cells throughout the entire mouse brain

Folder content:

- **figure_S03_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_S03_brainRenders.ipynb.
- **figure_S03_brainrenders.ipynb** heatmaps representing on brain coronal sections PV diffuse fluorescence and PV energy at mid-ontology resolution.
- **figure_S03_mainVisualizations.ipynb** all the other plots of the figure (scatterplots, heatmaps), representing the relationship of PV diffuse staining and PV energy as well as raw data of individual animals.

## Figure S04

### Colocalization of PNNs and PV cells in medium-resolution brain areas grouped by their major subdivision

Folder content:

- **figure_S04_colocalization.ipynb** Barplots showing colocalization of PNNs and PV cells at mid-ontology level.

## Figure S05

### Energy of PV<sup>-</sup> PNNs for coarse and medium-resolution brain areas

Folder content:

- **figure_S05_PVnegativePNNenergy.ipynb** Barplots showing PV<sup>-</sup> PNNs at coarse mid-ontology level.

## Figure S06

### WFA diffuse fluorescence and PV staining metrics in sensory cortical areas

Folder content:

#### WFA Diffuse Fluorescence in primary vs secondary areas by layers

- **figure_S06_WFAdiff_bylayer_.ipynb** Barplots showing WFA in sensory areas.

#### PV cell distribution in sensory cortical areas

- **figure_S06_primarySecondaryPV.ipynb** Barplots showing PV expression in sensory areas.

## Figure S07

### Determinants of PNN expression in the cortex

Folder content:

#### PV cell intensity and colocalization with PNNs in the sensory areas of the cortex

- **figure_S07_colocalizationPrimary_secondary.ipynb** Barplots showing PNN/PV colocalization in sensory areas.
- **figure_S07_PVPNN_primarySecondary.ipynb** Barplots showing the distribution of PV cells by intensity class in sensory areas.

#### Thalamic inputs from the association-cortex-related portion of the thalamus (DORpm) do not correlate with PNNs in sensory cortices

- **figure_S07_nonsensoryThalamus.ipynb** scatterplots representing the relationship between PNN expression in layers of the sensory cortex and density of DORpm afference.

#### Properties of PV cells in high-WFA and low-WFA cortical subnetworks

- **figure_S07_functionalGroups_PV_coloc.ipynb** Barplots showing PNN/PV colocalization and PV energy in cortical subnetworks.

## Supplementary data

Supplementary data are provided as .xlsx files produced with the following notebooks:

- **data_SD1_SD2.ipynb** whole-brain PNN/PV metrics.
- **data_SD3.ipynb** whole-brain PNN/PV colocalization metrics.
- **data_SD4.ipynb** correlation of staining metrics with gene expression.
