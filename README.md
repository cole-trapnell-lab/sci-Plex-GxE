# sci-Plex-GxE
sci-Plex-GxE is a platform for combined single-cell genetic and chemical screening at scale. This github repository documens processing pipelines and analysis scripts, including files used for data generation, analysis and interpretation.

### Contents

##### Data and data processing scripts
The data presented in our manuscript are derived from one of 4 single-cell screens. All raw and processed files can be found on the National Center for Biotechnology Information Gene Expression Ominibus (NCBI GEO) repository under series GSE225775.

| Name        | Experiment           |GEO Accession  | Analysis Scripts |
| :-------------: |:-----------:| :----:| :---:|
| sciPlexGxE_1      | HPRT1 & MMR  GxE Screens| [GSM7056148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056148)| HPRT1_MMR_GxE_screen|
| sciPlexGxE_2| Kinome GxE Screen      |  [GSM7056149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056149) | Kinome_GxE_screen|
| sciPlex_3 | GSC Chemical Genomic Screen      | [GSM7056150](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056150) | GSC_chemical_screen |
| sciPlex_4 | Combinatorial Chemical Genomic Screen      | [GSM7056151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056151)| combinatorial_chemical_screen |

##### Data processing and analysis
The scripts and code to process the data can be found in the [process\_from_raw](https://github.com/cole-trapnell-lab/sci-Plex-GxE/tree/main/process_from_raw) folder. The scripts and code to analyze the data can be found in the folders for each experiment. To access all files and functions used in this study, also download the [sci-Plex](https://github.com/cole-trapnell-lab/sci-plex) github repository from [Srivatsan, McFaline-Figueroa and Ramani et al. Science 2020](https://www.science.org/doi/10.1126/science.aax6234?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).
