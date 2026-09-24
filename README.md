# Trunk-River-Chlorobi

Script and data for: 

M.A. Moynihan, O.L. Mathieson, C.A. Crowley, E. Greene, H. Vanderscheuren, D. Dumit, R. Weed, K. Koop-Jakobsen, M. Kleiner, S.E. Ruff. Anoxygenic phototrophic Chlorobi use broad metabolic and resource acquisition strategies to support stable near-clonal blooms

Weed, R., Moynihan, M. A., Greene, E., Mathieson, O. L., Crowley, C.A., Junkins, E. N., Vanderscheuren, H., Jivaji, A. M., Kleiner, M., Ruff, S. E., Bloom-forming anoxygenic phototrophs heavily invest in anti-phage defense

## Moynihan et al Files

### Amplicon Scripts and Plots

### LoopSeq Scripts and Plots

### Metaproteomics 
[Metaproteomics](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Moynihan_et_al/Metaproteomics) folder contains scripts and files to make plots of metaproteomics-based biomass, community-wide sulfur and nitrogen cycle genes, and organism-normalized anvio plot of the proteome of the Prosthecochloris GSB-TRL01 metagenome-assembled genome over the sampling time series. 

## Weed et al Files
[Weed_et_al](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al): Contains code for Bloom-forming anoxygenic phototrophs heavily invest in anti-phage defense with the following sub-directories:

[viral_workflow](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al/viral_workflow): This folder contains the scripts for the generation of vMAGs, starting from raw reads. Although not designed to be run continuously, they were run in the following order:

* viral_assembly.sh
* viral_identification.sh
* viral_QC_and_combine.sh
* vMAG_generation.sh
* host_matching.sh

[bloomer_review_code](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al/bloomer_review_code):  This folder contains scripts for the review of bloom-forming organisms across the tree of life. Although not designed to be run continuously, they were run in the following order:
* download_and_filter_initial_sandpiper_db.sh
* sandpiper_db_filtering.R
* bloomer_and_rand_subset_defense_id.sh

[Figures_and_supplementary_tables](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al/Figures_and_supplementary_tables): This folder contains scripts for the production of figures 3, 5, and 6a, supplementary figures S2, S4, and S5-8, and supplementary tables 3, and 5-7.


#### additional files
[nanopore_sequence_processing](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al/nanopore_sequence_processing.sh): contains the workflow for processing the Nanopore data.

[read_mapping](https://github.com/moyn413/TrunkRiverChlorobi/tree/main/Weed_et_al/read_mapping.sh): contains code for mapping reads to MAGs and vMAGs for relative abundance estimations using CoverM.
