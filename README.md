# Human Lung Cell Atlas
![](https://s3.amazonaws.com/proddata.sagebase.org/3398210/3f086811-8854-43e7-b7b3-69d55cad900c/HLCA.jpg?response-content-disposition=attachment%3B%20filename%3D%22HLCA.jpg%22%3B%20filename%2A%3Dutf-8%27%27HLCA.jpg&response-content-type=image%2Fjpeg&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEFYaCXVzLWVhc3QtMSJHMEUCIQCZ%2BRVQl6W4Twc7fEKPmcSB%2B7C8%2Bs3vpJ4yZ9oGZvRNfgIgUmaI9IpmkJknYU5IsQ4T8okjnG2ST21OPhFS3AvR49wqtAMIHxAAGgwzMjU1NjU1ODU4MzkiDHyP4hYt726h2kj7YSqRA3k5GgPTcSdMLvssFp7gz48ia5doR3iuW6sdsVeKAN4LCKSktUwkwF8dvjOInVtCbHahoJO1w9zCuog9ObescUJ6eZgg4Eftlhsks9vgcqvjImdT5fV9uZajicvZYA15jM4GicDXPgw2zVpey40Fihd33Da%2FlRxjev6gLS5XYV7nkoq%2F0JZg5TzCUjoOS9CVwKPAFnEyeNFM0v59GAFGZaW5aPS%2BfQVnP2PDZTYuOc8hGTP%2FOEy1%2F4lyo09In%2FzFdnTHJKJ6CTqPJGYKGNJ4NoCd3K6ku2lMucCI2tGPjXvlCXhZEvaivLlVG4FYufLfSZtw4St%2FWlMpVf6SSKMfYQFg3l5fi7fPWDgApTwlUfashPQ7eaIYg3oMgKbWT0n4f8gqAuU0yN4%2FyJq4y3A57X96qLaXfn%2BqZaueTLiExWpg0cjp7yPKSNsQ9Vy1iExYZKLAqYWGKS4i24lL6YKoZ9g5C0ye7%2F84r1ON5R2F34AIl9lrdWMexRwz0dr0VofO7VeCEklAq2fc07fwHhuvn3uQMN34u%2FIFOusBRxGOhT0xCU5q%2F0adZm3i%2FgKc4Ctey10TfQ61DpO9Ez11cmE8YIvB%2F83sg6js8VFFOl3a3lPYHKqa9UKBWbwQlmMw%2FU0UBVgqA4M9vz1mn7B79IVkORh%2BMIoi1X%2Fx%2FWiGSrgf07BY%2FYntnKzCBENoon6i3gvxsCZq49e2Z7Uyon3PwRn44Wn4Wd12j3jBDIvGZxa7MdpnM0KXlYKEf2t1Rmg0zkOP6JYfIBOrs853C7TaXYEJPpsxwtlNvHcJHhceP59ZJB7og%2FPAq3aqtYEQLSHOp6T%2Bu44VwOi4d5Zs0ff5aE978YBn3KcxwA%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200220T230554Z&X-Amz-SignedHeaders=host&X-Amz-Expires=29&X-Amz-Credential=ASIAUXTJYTGX6ZJBB7NM%2F20200220%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=6a0ab7083766964445aa1cc9392c1a83bb1ea8c375978bae55ddfbb20e60227a)
Although single cell RNA sequencing studies have begun providing compendia of cell expression profiles, it has proven more difficult to systematically identify and localize all molecular cell types in individual organs to create a full molecular cell atlas. From droplet- and plate-based single cell RNA sequencing applied to ~75,000 human lung and blood cells, combined with a multi-pronged cell annotation approach that includes extensive tissue localization, we have defined the gene expression profiles and anatomical locations of 58 cell populations in the human lung, including 41 of 45 previously known cell types or subtypes and 14 new ones. Learn more in our [manuscript](https://www.biorxiv.org/content/10.1101/742320v1) or explore the data in your browser with [cellxgene](http://hlca.ds.czbiohub.org).

### Annotation
These R markdown notebooks import patient-specific gene count/UMI tables (SS2/10x) with [Seurat](https://satijalab.org/seurat/), seperate the cells by tissue compartment, iteratively subclsuter them, and then assign each biologically meaningful cluster an identity based on canonical marker genes, novel bulk RNAseq markers (immune cells), tissue location, and biochemical function. Caution: The specific version of Seurat used (as well as its dependencies) can produce slightly different clustering results. Early steps in these notebooks remove cells from patient-matched diseased regions of the lung, which are not yet released.

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
