# wholeBrain_PNN_analysis
Analysis code to reproduce all figures of the paper _A Comprehensive Atlas of Perineuronal Net Distribution and Colocalization with Parvalbumin in the Adult Mouse Brain_, published in [journal name]() <or in bioRxiv>


### Dependencies

To run the scripts in this repository you will need two conda environments.
- An environment with the packages listed in requirements.txt (python version: 3.6 or 3.7)
- An environment with the packages listed in requirements_brainrender.txt

The second environment is needed only for the renderings of brain slices, realised through [brainrender](https://github.com/brainglobe/brainrender) and [bg-heatmaps](https://github.com/brainglobe/bg-heatmaps)

Thus, unless otherwise stated, run the code in the first environment.

## **wholebrain_tools package**
---
Analysis in this work were conducted with the help of custom code wrapped in a python package. This package contains the following modules:
- **aba**: handling of the Allen Brain Atlas 
- **stats**: tools for statistical analysis
- **graphics**: tools for custom plots
- **dataIO**: file handling

Importantly, the repository contains a path.ini.

## **Figure 2**

### Distribution of WFA-positive PNNs throughout the entire mouse brain
---
Plots in this figures provide a graphical visualization of the distriburion of WFA<sup>+</sup> perineuronal nets throughout the mouse brain.

The folder contains the following notebooks:
- **figure_S03_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_02_brainRenders.ipynb
- **figure_02_brainRenders.ipynb** heatmaps representing on brain coronal sections WFA diffuse fluorescence and PNN energy at mid-ontology resolution. Note: this notebook needs a conda environment with the packages listed in requirements_brainrender.txt.
- **figure_02_mainVisualizations.ipynb** all the other plots of the figure (scatterplots, heatmaps), representing the relationship of WFA diffuse staining and PNN energy as well as raw data of individual animals.

## **Figure 3**

### Brain-wide interactions between PNNs and PV cells
---
Plots in this figures provide a graphical visualization of the colocalization of PNNs and PV cells throughout the mouse brain.

The folder contains the following notebooks:
- **figure_03_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_02_brainRenders.ipynb
- **figure_03_colocalizationBrainRenders.ipynb** heatmaps representing on brain coronal sections the percentage of PNNs surrounding a PV cells and the percentage of PV cells ensheated by PNNs. Note: this notebook needs a conda environment with the packages listed in requirements_brainrender.txt.
- **figure_03_Correlation_pv_wfa.ipynb** scatterplots representing the relationship of WFA diffuse staining/PNN energy and PV energy.
- **figure_03_colocalization.ipynb** barplots to visulize PNN/PV colocalization at the coarse level

## **Figure 4**

### PNNs aggregation depends on PV expression levels
---
Plots in this figures explore the relationship between probability of having a PNN and PV expression levels throughout the mouse brain.

The folder contains the following notebooks:
- **figure_04_colocProbability.ipynb**  scatterplots representing the relationship between PNN expression and intensity staining of individual PV cells.

## **Figure 5**

### Principles of PNN organization in cortical areas
---
Plots in this figures explore the expression of PNNs in the cortex.

The folder contains the following notebooks:
- **figure_05_connectome.ipynb** scatterplots representing the relationship between PNN expression in layers of the sensory cortex and density of thalamic afference.
- **figure_05_corticalHeatmap.ipynb** heatmaps representing PNN expression in cortical areas.
- **figure_05_functionalRegions.ipynb** scatterplots representing the relationship of PNN expression and intensity staining of individual PV cells.
- **figure_05_primarySecondary.ipynb** barplots representing PNN expression in primary vs associative sensory areas.

## **Figure 6**

### Gene expression correlates of PNNs
---
Plots in this figures explore correlation of PNN energy with known PNN molecular markers and the biological processes in which the genes correlated or anticorrelated with WFA staining are involved.

- **figure_06_downloadGeneExprABA.ipynb** download of the AGEA dataset as a .csv file (see [Lein et al., 2007](https://doi.org/10.1038/nature05453)).
- **figure_06_geneOntology.ipynb** production of gene lists for gene ontology analysis and plotting of the results.
- **figure_06_matrisome.ipynb** production of gene lists for overrepresentation analysis in the matrisome dataset [Naba et al., 2016](https://doi.org/10.1016/j.matbio.2015.06.003).
- **figure_06_plotMarkers.ipynb** scatterplots representing PNN/PV levels as a function of the expresssion of known molecular markers.
- **figure_06_preprocess_dataWFA.ipynb** correlation analysis of WFA fluorescence/PNN energy and gene expression data of AGEA dataset.
- **figure_06_preprocess_dataPV.ipynb** correlation analysis of PV energy and AGEA dataset.

## **Figure S01**

### The scores assigned by the scoring models correlate with raters’ agreement
---
Plots in this figures show the performances of our scoring model

The folder contains the following notebooks:
- **figure_S01_networkMetrics.ipynb** boxplots representing the score assigned by the model to cells at different levels of agreement

## **Figure S02**

### PNN energy and WFA diffuse fluorescence measurements for medium-resolution brain areas grouped by their major subdivision
---
Plots in this figures show PNN expression metrics at mid-ontology resolution

The folder contains the following notebooks:
- **figure_S02.ipynb** Barplots showing PNN energy and WFA diffuse fluorescence at mid-ontology level.

## **Figure S03**

### Distribution of PV-positive cells throughout the entire mouse brain
---

The folder contains the following notebooks:
- **figure_S03_prepareDataForBrainRender.ipynb** preprocessing of data required by figure_S03_brainRenders.ipynb
- **figure_S03_brainrenders.ipynb** heatmaps representing on brain coronal sections PV diffuse fluorescence and PV energy at mid-ontology resolution. Note: this notebook needs a conda environment with the packages listed in requirements_brainrender.txt.
- **figure_02_mainVisualizations.ipynb** all the other plots of the figure (scatterplots, heatmaps), representing the relationship of PV diffuse staining and PV energy as well as raw data of individual animals.

## **Figure S04**

### Colocalization of PNNs and PV cells in medium-resolution brain areas grouped by their major subdivision
---
Plots in this figures show PNN/PV colocalization metrics at mid-ontology resolution

The folder contains the following notebooks:
- **figure_S04_colocalization.ipynb** Barplots showing colocalization of PNNs and PV cells at mid-ontology level.

## **Figure S05**

### WFA Diffuse Fluorescence in primary vs secondary areas by layers
---
Plots in this figures show WFA diffuse fluorescence in sensory primary versus associative areas, split by layers.

The folder contains the following notebooks:
- **figure_S05_WFAdiff_bylayer_.ipynb** Barplots showing WFA in sensory areas.

## **Figure S06**

### PV cell distribution in sensory cortical areas.
---
Plots in this figures show PV distribution the sensory areas of the cortex.

The folder contains the following notebooks:
- **figure_S06_colocalizationPrimary_secondary.ipynb** Barplots showing PNN/PV colocalization in sensory areas.
- **figure_S06_primarySecondaryPV.ipynb** Barplots showing PV expression in sensory areas.

## **Figure S08**

### Properties of PV cells in high-WFA and low-WFA cortical subnetworks
---
Plots in this figures show properties of PV cells in cortical subnetworks

The folder contains the following notebooks:
- **figure_S08_functionalGroups_PV_coloc.ipynb** Barplots showing PNN/PV colocalization and PV energy in cortical subnetworks.


## **Supplementary data**
---
Supplementary data are provided as .xlsx files produced with the following notebooks:

- **data_SD1_SD2** whole-brain PNN/PV metrics
- **data_SD3** whole-brain PNN/PV colocalization metrics
- **data_SD4** correlation of staining metrics with gene expression