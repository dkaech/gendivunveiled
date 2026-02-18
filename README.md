This repository contains data and code used in the following paper:

**Genetic diversity unveiled: cost-effective methods for grassland species** 

Damian Käch, Miguel Loera-Sánchez, Bruno Studer,  Roland Kölliker

Preprint 2025, doi: https://doi.org/10.1101/2025.07.23.666374

Peer-reviewed article 2026, doi: https://doi.org/10.1111/1755-0998.70108


## 🧪 Overview and Materials and Methods
Permanent grasslands are the basis for sustainable ruminant livestock production and provide various ecosystem services. They are mainly composed of outcrossing plant species, leading to populations with high genetic diversity (i.e., intraspecific diversity). Grasslands of high plant genetic diversity (PGD) can better cope with environmental stress and have stabilised biomass productivity. Additionally, they are valuable reservoirs of genetic resources used for forage plant breeding. To detect undesired changes and intervene accordingly, monitoring PGD in these grasslands is key. Despite the availability of various molecular genetic approaches, PGD monitoring is often neglected in biodiversity reports, which is attributed to a lack of standardised and affordable indicators of genetic diversity in natural populations.

To assess PGD of agronomically relevant grassland species, we applied multispecies amplicon sequencing (MSAS) and genotyping-by-sequencing (GBS), resulting in three data sets. Using MSAS, we analysed 39 samples based on five species (*Dactylis glomerata* L., *Festuca pratensis* Huds., *Lolium perenne* L., *Trifolium pratense* L., *Trifolium repens* L.). The sample set contains 30 single-accession (SA) seedling samples (five species, two accessions (A and B) per species, three replicates) and nine mixed-species (MS) seedling samples (three compositions, three replicates). The latter were prepared by pooling DNA from SA samples for three different compositions: 1) MS-A100, containing accession A of each species at equal amounts, 2) MS-B100, containing accession B of each species at equal amounts, and 3) MS-AB50, for equal amounts of MS-A100 and MS-B100. 

Furthermore, we prepared an extended *L. perenne* sample set consisting of 42 samples based on six cultivars. This sample set contained 18 single-cultivar samples (six cultivars, three replicates) and 18 mixtures of two cultivars (three mixtures, two mixing ratios (50:50 and 75:25), three replicates). In addition to these samples, which were based on a greenhouse experiment, the sample set contained six 50:50-ratio mixtures based on a field experiement (one mixture, two locations, three replicates). Subsequently, the sample set was analysed using MSAS and GBS.

## 📁 Repository Structure
```
├── data/
│   ├── ADMIXTURE/
│   │   ├── GBS/
│   │   │   ├── bootstrapping/
│   │   │   ├── focus_groups/
│   │   │   └── gbs_topweide_nolt.list
│   │   ├── MSAS/
│   │   │   ├── bootstrapping/
│   │   │   └── focus_groups/
│   ├── BAM_headers/
│   |   ├── GBS/
│   |   └── MSAS/
│   ├── MSAS reference sequences/
│   ├── STRUCTURE/
│   │   ├── GBS/
│   │   └── MSAS/
│   ├── VCF/
│   │   ├── GBS/
│   │   └── MSAS/
├── results/
│   ├── figures/
│   └── tables/
├── scripts/
│   ├── 01_poolfstat.R
│   ├── 02_poolfstat_permanova.R
│   ├── 03_taxsort_poolfstat.R
│   ├── 04_fst_pairwise_poolfstats.R
│   ├── 05_taxsort_poolfstat_DAPC.R
│   ├── 06_taxsort_pure_and_mixed_Fst.R
│   ├── 07_taxonomic_assignment_rate.R
│   ├── 11_init.R
│   ├── 12_GBS_loadData_FST_AF.R
│   ├── 13_GBS_DAPC.R
│   ├── 14_GBS_ADMIXTURE_STRUCTURE.R
│   ├── 15_MSAS_loadData_FST_AF.R
│   ├── 16_MSAS_DAPC.R
│   ├── 17_MSAS_ADMIXTURE_STRUCTURE.R
│   └── 18_GBSandMSAS.R
```


## 📊 Data Overview

- `ADMIXTURE`: Files for ADMIXTURE analysis of GBS and MSAS data based on extended *L. perenne* sample set. `gbs_topweide_nolt.list` contains the sample names, `focus_groups` contains files with cluster memberships, and `bootstrapping` contains files with standard errors for cluster memeberships.
- `STRUCTURE`: Files for STRUCTURE analysis of GBS and MSAS data based on extended *L. perenne* sample set. Files with the suffix `f` contain cluster memberships of pure and `gbs_topweide_priors_man.txt` of the mixed samples.
- `VCF`: VCF files of the multispecies (`MSAS/MS/`, `MSAS/SA/`, `MSAS/SAMS/`) and the *L. perenne* (`GBS/`, `MSAS/extended_loper/`) sample set.
- `BAM_headers`: Headers of the BAM files of the multispecies (`MSAS/MS/`, `MSAS/SA/`, `MSAS/SAMS/`) and the *L. perenne* (`GBS/`, `MSAS/extended_loper/`) sample set.
- `MSAS reference sequences`: Reference sequences used for multispecies sample set.


## 📈 Results Overview
Figures and tables/dataframes resulting from the scripts separated by sample set type.

## 📜 Scripts Overview

Scripts with the prefix `0` were used for the multispecies samples set and for computing fixation index values (*F*<sub>ST</sub>) PERMANOVA of all sample sets. Scripts with the prefix `1` were used for the extended *L. perenne* sample set. Scripts with the same prefix have to be run sequentially (i.e., in order of increasing second number in file name).

- `01_poolfstat.R`: Summaryse poolfstat results for multispecies samples.
- `02_poolfstat_permanova.R`: Analyse multispecies samples using poolfstat and PERMANOVA.
- `03_taxsort_poolfstat.R`: Analyse taxonomically sorted reads.
- `04_fst_pairwise_poolfstats.R`: Analyse correlation of poolfstat results between MSAS and GBS (based on extended *L. perenne* sample set).
- `05_taxsort_poolfstat_DAPC.R`: Analyse allele frequencies using discriminant analysis of principal components (DAPC; based on multispecies data set).
- `06_taxsort_pure_and_mixed_Fst.R`: Analyse SAMS samples of multispecies sample set using poolfstat and PERMANOVA).
- `07_taxonomic_assignment_rate.R`: Estimate taxonomic assignment rate.
- `11_init.R`: Initialise by defining paths, loading packages and functions (for analysis of extended *L. perenne* sample set).
- `12_GBS_loadData_FST_AF.R`: Load VCF file and compute *F*<sub>ST</sub> and allele frequencies (for analysis of GBS data based on extended *L. perenne* sample set).
- `13_GBS_DAPC.R`: Analyse allele frequencies using a discriminant analysis of principal components (DAPC; for analysis of GBS data based on extended *L. perenne* sample set).
- `14_GBS_ADMIXTURE_STRUCTURE.R`: Analyse population structure using ADMIXTURE and STRUCTURE (for analysis of GBS data based on extended *L. perenne* sample set).
- `15_MSAS_loadData_FST_AF.R`: Load VCF file and compute *F*<sub>ST</sub> and allele frequencies (for analysis of MSAS data based on extended *L. perenne* sample set).
- `16_MSAS_DAPC.R`: Analyse allele frequencies using DAPC (for analysis of MSAS data based on extended *L. perenne* sample set).
- `17_MSAS_ADMIXTURE_STRUCTURE.R`: Analyse population structure using ADMIXTURE and STRUCTURE (for analysis of MSAS data based on extended *L. perenne* sample set).
- `18_GBSandMSAS.R`: Compare GBS and MSAS data based on extended *L. perenne* sample set.

