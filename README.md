# High-throughput RNA isoform sequencing using programmable cDNA concatenation

## Abstract
Alternative splicing is a core biological process that enables profound and essential diversification of gene function. Short-read RNA sequencing approaches fail to resolve RNA isoforms and therefore primarily enable gene expression measurements - an isoform unaware representation of the transcriptome. Conversely, full-length RNA sequencing using long-read technologies are able to capture complete transcript isoforms, but their utility is deeply constrained due to throughput limitations.  Here, we introduce _**MAS-ISO-seq**_, a technique for programmably concatenating cDNAs into single molecules optimal for long-read sequencing, _**boosting the throughput >15 fold to nearly 40 million cDNA reads per run on the Sequel IIe sequencer**_. We validated unambiguous isoform assignment with MAS-ISO-seq using a synthetic RNA isoform library and applied this approach to single-cell RNA sequencing of tumor-infiltrating T cells. Results demonstrated a >30 fold boosted discovery of differentially spliced genes and robust cell clustering, as well as canonical PTPRC splicing patterns across T cell subpopulations and the concerted expression of the associated hnRNPLL splicing factor. Methods such as MAS-ISO-seq will drive discovery of novel isoforms and the transition from gene expression to transcript isoform expression analyses.

## Authors
_**Aziz M. Al’Khafaji**_<sup>1*†</sup>, _**Jonathan T. Smith**_<sup>1*</sup>, _**Kiran V Garimella**_<sup>1*†</sup>, _**Mehrtash Babadi**_<sup>1*†</sup>, Moshe Sade-Feldman<sup>1,2</sup>, Michael Gatzen<sup>1</sup>, Siranush Sarkizova<sup>1</sup>, Marc A. Schwartz<sup>1,3,4</sup>, Victoria Popic<sup>1</sup>, Emily M. Blaum<sup>1,2</sup>, Allyson Day<sup>1</sup>, Maura Costello<sup>1</sup>, Tera Bowers<sup>1</sup>, Stacey Gabriel<sup>1</sup>, Eric Banks<sup>1</sup>, Anthony A. Philippakis<sup>1</sup>, Genevieve M. Boland<sup>5</sup>, Paul C. Blainey<sup>1,6,8,†</sup>, Nir Hacohen<sup>1,7,10,11,†</sup>

1. Broad Institute of Harvard and MIT, Cambridge, MA, USA
2. Department of Medicine, Center for Cancer Research, Massachusetts General Hospital, Boston, MA, USA
3. Department of Pediatrics, Harvard Medical School, Boston, Massachusetts, USA.
4. Division of Hematology/Oncology, Boston Children's Hospital, Boston, Massachusetts, USA.
5. Department of Pediatric Oncology, Dana Farber Cancer Institute, Boston, Massachusetts, USA.
6. Division of Surgical Oncology, Massachusetts General Hospital, Harvard Medical School, Boston, MA
7. Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA
8. Center for Cancer Research, Massachusetts General Hospital and Harvard Medical School, Boston, MA, USA
9. Koch Institute for Integrative Cancer Research at Massachusetts Institute of Technology, Cambridge, MA, USA
10. Harvard Medical School, Boston, MA, USA
11. Center for Immunology and Inflammatory Diseases, Massachusetts General Hospital, Charlestown, MA, USA

\* - These authors contributed equally  
† - Corresponding authors


## Data

- All data from this study are available online (or are in the process of being uploaded).  

There were two datasets from this study: 

| Dataset | Number of Samples | Location |
|---|---|---|
| Human tumor-infiltrating CD8+ T cells | 2 | _Release in progress..._ |
| Spike-in RNA Variant Control Mix data ([SIRVs set 4](https://www.lexogen.com/sirvs/sirv-sets/), [Lexogen](https://www.lexogen.com/)) | 2 | [Terra](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/MAS-seq%20-%20Data%20Segmentation%20and%20Alignment/data)* |

\* - The SIRV samples were prepared with two library preparation techniques: a length 10 MAS-ISO-seq array and a length 15 MAS-ISO-seq array.  They were multiplexed into a single pooled sample and sequenced in a single run on a PacBio Sequel IIe.  Our software package, [Longbow](https://github.com/broadinstitute/longbow/releases/tag/v0.2.2), was then used to demultiplex the single SIRV multiplexed sample into two outputs - one for the length 15 array and one for the length 10 array.  These demultiplexed files are what is currently available in the Terra workspace.

## Terra Workspace Example
A [Terra](https://terra.bio) workspace with an example of how to process MAS-ISO-seq data can be found here:

- [MAS-seq - Data Segmentation and Alignment](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/MAS-seq%20-%20Data%20Segmentation%20and%20Alignment)

This workspace is an example of how to segment and align MAS-ISO-seq data.  

The data in this workspace are the same Spike-in RNA Variant Control Mix ([SIRVs set 4](https://www.lexogen.com/sirvs/sirv-sets/), [Lexogen](https://www.lexogen.com/)) samples that we used as controls in the paper.

## Pre-print of the Paper
A preprint of the paper can be found on bioRxiv here: _LINK TO BE ADDED_

