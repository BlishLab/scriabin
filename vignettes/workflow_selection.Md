

## Choosing the correct Scriabin workflow ##

The graphic below assists users in selecting an appropriate workflow for Scriabin

![](vignettes/workflow_selection.png)

The first step in selecting a workflow is considering the magnitude of the dataset being analyzed -- or more, specifically, what number of cell-cell pairs would be generated during the CCIM workflow? Because the CCIM scales exponentially with dataset size, it is frequently impractical to generate a CCIM for an entire dataset. For instance, a sample with 1,000 cells would generate a CCIM with 1 million observations. Determining what dataset size is too big for the CCIM depends both on the magnitude of interaction between each pair of cells, as well as the computational resources available to the user. We generally find that CCIM with a total cell-cell pair count of <500,000 are easily amenable to analysis locally. 

Note: we recommend using Scriabin's CCIM workflow to analyze pairs of that were actually capable of interacting in reality. For example, consider a PBMC dataset of 2k cells, where 1k cells each come from two different donors. We only generate a CCIM *per sample*, so the CCIM size for this dataset would be 1 million for each donor, or 2 million cell-cell pairs in total. 
