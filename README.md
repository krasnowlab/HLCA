# Human Lung Cell Atlas
![](https://hlca.ds.czbiohub.org/images/HLCA.jpg)
Although single cell RNA sequencing studies have begun providing compendia of cell expression profiles, it has proven more difficult to systematically identify and localize all molecular cell types in individual organs to create a full molecular cell atlas. From droplet- and plate-based single cell RNA sequencing applied to ~75,000 human lung and blood cells, combined with a multi-pronged cell annotation approach that includes extensive tissue localization, we have defined the gene expression profiles and anatomical locations of 58 cell populations in the human lung, including 41 of 45 previously known cell types or subtypes and 14 new ones. Learn more in our [manuscript](https://www.biorxiv.org/content/10.1101/742320v1) or explore the data in your browser with [cellxgene](http://hlca.ds.czbiohub.org).

### Annotation
These R markdown notebooks import patient-specific gene count/UMI tables (SS2/10x) with [Seurat](https://satijalab.org/seurat/), seperate the cells by tissue compartment, iteratively subclsuter them, and then assign each biologically meaningful cluster an identity based on canonical marker genes, novel bulk RNAseq markers (immune cells), tissue location, and biochemical function. Note: The specific version of Seurat used (as well as its dependencies) can produce slightly different clustering results. Early steps in these notebooks remove cells from patient-matched diseased regions of the lung, which are not yet released.

### Analysis
These notebooks import annotated, patient-specific Seurat objects (see Data availability below) from the annotation notebooks and merges them for downstream analyses. They explore the biochemical functions of lung cell types and the cell-selective transcription factors and optimal markers for making and monitoring them; define the cell targets of circulating hormones and predicts local signaling interactions including sources and targets of chemokines in immune cell trafficking and expression changes on lung homing; and identify the cell types directly affected by lung disease genes. They also compare human and mouse cell types to identify cell types whose expression profiles have been substantially altered by evolution, revealing extensive plasticity of cell-type-specific gene expression. These notebooks output the Figures and Tables used in the manuscript.

There are several files used by the Analysis notebooks that are too big for GitHub. To obtain them:

(1) Clone the GitHub repository

```
git clone http://github.com/krasnowlab/hlca
```

(2) Setup the Synapse [CLI](https://python-docs.synapse.org/build/html/CommandLineClient.html)

(3) Navigate to the repository

(4) Run the download script

```
bash fetchData.bash
```


### Barga et al and Reyfman et al
We reannotated two existing human lung single cell RNAseq datasets from the [Teichmann](https://www.nature.com/articles/s41591-019-0468-5) and [Misharin](https://www.atsjournals.org/doi/full/10.1164/rccm.201712-2410OC) labs using our approach. Notebooks and CSVs containing the new annotations are provided.

### Data availability
Download count/UMI tables, metadata, and Seurat and scanpy objects for your own bioinformatic piplines from [Synapse](https://www.synapse.org/#!Synapse:syn21041850/wiki/600865). Sequencing reads will be deposited on the European Genome-phenome Archive following final approval of our Data Access Agreement and Committee by Stanford's Privacy Office.
