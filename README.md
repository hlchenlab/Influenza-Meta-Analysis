# Meta-analysis of Influenza studies

## Introduction:
The Influenza virus is one of the most pathogenic viruses in the world and has previously shown the potential to establish worldwide pandemics throughout its troubled history. Much of the current literature has described in much detail the intricacies of how influenza infection can interact with cellular processes for the benefit of viral replication. One such interaction is the well-studied subversion of the host immune response to the invading influenza virus in infected cells. Several mechanisms have been proposed to be essential for this trait, but many have proved to be very dependent to the conditions used in individual infection experiments. Furthermore, the genes identified as being targeted by influenza specific signalling to interfer with intracellular activation of antiviral activity has shown much inconsistancy between individual studies. We describe in this repository a unique methodology to robustly identify recurrent gene signatures across multiple influenza studies thereby removing study specific gene expression and enhancing the expression signal induced by influenza infection. 

The complexity of genomes across many multi-cellular organisms is often highlighted by how genes function together in highly controlled gene networks. Therefore, components of infected cells, there are identified as targets of influenza interference will ultimately result in an altered expression profile of select genes that will correspond to altered functional activity of these gene networks. Coexpression analysis can be used to identify how genes correlate in their expression profiles and identify gene network modules within diffential expression gene lists. We utilise a two-pronged approach in characterizing gene expression as a result of influenza specific-signalling processes by using meta-analysis and co-expression. The end-result is the identification gene co-expression modules that represent universal functional gene co-expression networks induced directly by influenza infection.

## Contents of this respository

### Aquisition of studies

RNA-seq studies were acquired from an extensive search on the Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo) database (Search for human studies was undertaken in March 2019 and search for mouse studies was done in November-December 2019). Studies were required to have at least two conditions of infection (i.e. Infected samples with a wildtype influenza viral strain and uninfected samples). In total, 8 human studies and 7 mouse were acquired. 

#### Human Studies

|GSE ID:|SRA ID:|Molecule:|Library Preparation:|Sample Number:|Citation:|
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
|GSE84204|SRP078169|Total RNA|Paired-end|15|Alculumbre et al., 2018|
|GSE89008|SRP091886|Total RNA|Paired-end|52*|Heinz et al., 2018|
|GSE97672|SRP103821|Total RNA|Paired-end|64*|Heinz et al., 2018|
|GSE103477|SRP117084|Total RNA|Paired-end/Single-end|20**|Heinz et al., 2018|
|GSE103604|SRP117055|Total RNA|Paired-end/Single-end|24*|Zhao et al., 2018|
|GSE104168|SRP118721|Total RNA|Single-end|30*|Forst et al., 2017|
|GSE156060|SRP277089|Total RNA|Paired-end|24**|This study|
|GSE156152|SRP277269|Total RNA|Paired-end|12*|This study|

\* Not all samples were used, samples derived from knockdown, transfection experiments/mutant infections/infection experiments under 3 hpi were excluded.

** Not all samples were used, samples derived from ChIP-seq and mutant virus infections were excluded

#### Mouse Studies

|GSE ID:|SRA ID:|Molecule:|Library Preparation:|Sample Number:|Citation:|
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
|GSE117029|SRP153200|Total RNA|Paired-end|25|Sengupta et al., 2019|
|GSE100522|SRP110570|Total RNA|Paired-end|9|Hu et al., 2018|
|GSE107488|SRP125820|poly A|Paired-end|14*|Yildiz et al., 2018|
||SRP061303|poly A|Paired-end|21*||
||ERP020504|poly A|Single-end|24*|Steed et al., 2017|
|GSE49933|SRP028868|poly A|Single-end|23|Altboum et al., 2014|
|GSE52405|SRP033021|Total RNA|Paired-end|78|Josset et al., 2014|

\* Not all samples were used, samples derived from non-lung tissue/transgenic mouse experiments/infection experiments under 1 dpi were excluded.

### Alignments

All alignments were performed using STAR (v2.6.0c) (Dobin et al., 2013) with default settings to acquire counts using either the human or mouse genomes. THe human genome and annotation were obtained from GENCODE (v29) (https://www.gencodegenes.org/human/release_29.html). Mouse genome and annotation were also acquired from GENCODE (v24) (https://www.gencodegenes.org/mouse/release_M24.html).

### Differential Expression

### Co-expression and Clustering Analysis

### Cross-species conservation

###

## References

