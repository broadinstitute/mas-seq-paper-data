# MAS-seq / Longbow Workflows

This folder contains the workflow for processing MAS-seq data with
Longbow and all other required tools to produce a transcript count
matrix for further analysis.

These [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) files have several numbered tasks
which are executed in approximate numbered order when submitted to the cloud via [Cromwell](https://cromwell.readthedocs.io/en/stable/) or [Terra](https://terra.bio/)
(the order is not exactly numerical because cromwell performs call tree analysis to 
parallelize tasks that can be run concurrently). 

## User Caveat
The workflows defined here should be considered an alpha/beta release.
They produce valid data, however they are not optimized for general use and should not be run without inspection and optimization. 

## Docker Images
All docker images used in these WDL files are stored in the 
[Google Container Registry](https://cloud.google.com/container-registry) for the DSP Methods Long Reads Methods and 
Applications (LRMA) group at The Broad Institute of MIT and Harvard.  The full paths are specified inline in each WDL task.

## Workflow Versions
The workflows in this paper were modified in the time between the preprint release and the final release.  They have both been preserved here. 

### Final Submission Workflows
The workflows in the [final](final) folder correspond to processing performed on the human peripheral blood mononuclear cell (PBMC) sample data in the revised submission of the paper.

NOTE: The SIRV analysis did not change from the preprint, so these WDLs are not reproduced in this folder.
NOTE: The downsampled WDL files are a subset of the Main MAS-ISO-seq analysis, so the inputs for those WDLs are not included here.

### Preprint Workflows
The workflows in the [preprint](preprint) folder correspond to processing performed on the human peripheral blood mononuclear cell (PBMC) sample data in the preprint submission of the paper.
