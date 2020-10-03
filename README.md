# Meta-analysis of Influenza studies

## Introduction:
The Influenza virus is one of the most pathogenic viruses in the world and has previously shown the potential to establish worldwide pandemics throughout its troubled history. Much of the current literature has described in much detail the intricacies of how influenza infection can interact with cellular processes for the benefit of viral replication. One such interaction is the well-studied subversion of the host immune response to the invading influenza virus in infected cells. Several mechanisms have been proposed to be essential for this trait, but many have proved to be very dependent to the conditions used in individual infection experiments. Furthermore, the genes identified as being targeted by influenza specific signalling to interfer with intracellular activation of antiviral activity has shown much inconsistancy between individual studies. We describe in this repository a unique methodology to robustly identify recurrent gene signatures across multiple influenza studies thereby removing study specific gene expression and enhancing the expression signal induced by influenza infection. 

The complexity of genomes across many multi-cellular organisms is often highlighted by how genes function together in highly controlled gene networks. Therefore, components of infected cells, there are identified as targets of influenza interference will ultimately result in an altered expression profile of select genes that will correspond to altered functional activity of these gene networks. Coexpression analysis can be used to identify how genes correlate in their expression profiles and identify gene network modules within diffential expression gene lists. We utilise a two-pronged approach in characterizing gene expression as a result of influenza specific-signalling processes by using meta-analysis and co-expression. The end-result is the identification gene co-expression modules that represent universal functional gene co-expression networks induced directly by influenza infection.

## Description of contents within this respository

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
|GSE49933|SRP028868|poly A|Single-end|23|Altboum et al., 2014|
|GSE52405|SRP033021|Total RNA|Paired-end|78|Josset et al., 2014|
|GSE100522|SRP110570|Total RNA|Paired-end|9|Hu et al., 2018|
|GSE107488|SRP125820|poly A|Paired-end|14*|Yildiz et al., 2018|
|GSE117029|SRP153200|Total RNA|Paired-end|25|Sengupta et al., 2019|
||ERP020504|poly A|Single-end|24*|Steed et al., 2017|
||SRP061303|poly A|Paired-end|21*||

\* Not all samples were used, samples derived from non-lung tissue/transgenic mouse experiments/infection experiments under 1 dpi were excluded.

### Alignment

All alignments were performed using STAR (v2.6.0c) (Dobin et al., 2013) with default settings to acquire counts using either the human or mouse genomes. THe human genome and annotation were obtained from GENCODE (v29) (https://www.gencodegenes.org/human/release_29.html). Mouse genome and annotation were also acquired from GENCODE (v24) (https://www.gencodegenes.org/mouse/release_M24.html). Bash script used for alignments are included.

### Differential Expression and meta-analysis

For each study, samples were classed into two categories, infected (samples generated from wildtype infections only), or mock infected (uninfected control). All other samples from studies were excluded. For human studies, genes were classified as differentially expressed (DE) if DESeq2 (Love et al., 2014) identified genes with a log2FC >= 1.5 and a p-adjusted value <= 0.05. Mouse studies used less stringent DE thresholds of log2FC >= 1 and p-adjusted <= 0.05. Rscript used for perform differential expression with DESeq2 and the generation of volcano-plots and heatmaps is included. For meta-analysis, genes identified as DE were given a score of 1 ortherwise given a 0 per study. Significant recurrence thresholds of gene differential expression was then calculated across all studies for each species. Differentially expressed genes (DEGs) with a recurrences at or above the calculated threshold where retained for clustering. Rscript for recurrence analysis is included in this repository. Accessory functions for this are available at https://github.com/sarbal/OutDeCo/tree/master/R. Human and mouse gene lists are also in this repository.

### Co-expression and Clustering Analysis

Rscript (Steps 1-4) describes construction of gene-gene co-expression network. In brief, studies had to have at least 10 samples to be included. Gene with significant expression (at least 10 reads, in at least 10 samples for 10 studies) were only included in the final network. Clustering_function.R is the function neccessary to perform hierarchical clustering on genes identified as statistically recurrent from meta-analysis. Study IDs used in human network construction are included in coexpression studies.csv. Construction of networks followed the basis of works done previously (Lee et al., 2020) by our collaborators at CSHL.

### Gene Enrichment Analysis

Code and annotation for Gene Ontology (http://geneontology.org/docs/downloads/, downloaded January 2020) is avaiable in this directory.

### Cross-species conservation

Conservation analysis used high confidence orthologues present in both human and mouse. Orthologue ids were downloaded from Ensembl (v100) (http://jan2020.archive.ensembl.org/index.html). An Rdata file is available here. Rscript for riverplots and Venn-diagram plots is also here. 

## References

ALCULUMBRE, S. G., SAINT-ANDRE, V., DI DOMIZIO, J., VARGAS, P., SIRVEN, P., BOST, P., MAURIN, M., MAIURI, P., WERY, M., ROMAN, M. S., SAVEY, L., TOUZOT, M., TERRIER, B., SAADOUN, D., CONRAD, C., GILLIET, M., MORILLON, A. & SOUMELIS, V. 2018. Diversification of human plasmacytoid predendritic cells in response to a single stimulus. Nat Immunol, 19, 63-75.

ALTBOUM, Z., STEUERMAN, Y., DAVID, E., BARNETT-ITZHAKI, Z., VALADARSKY, L., KEREN-SHAUL, H., MENINGHER, T., MENDELSON, E., MANDELBOIM, M., GAT-VIKS, I. & AMIT, I. 2014. Digital cell quantification identifies global immune cell dynamics during influenza infection. Mol Syst Biol, 10, 720.

DOBIN, A., DAVIS, C. A., SCHLESINGER, F., DRENKOW, J., ZALESKI, C., JHA, S., BATUT, P., CHAISSON, M. & GINGERAS, T. R. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29, 15-21.

FORST, C. V., ZHOU, B., WANG, M., CHOU, T. W., MASON, G., SONG, W. M., SCHADT, E., GHEDIN, E. & ZHANG, B. 2017. Integrative gene network analysis identifies key signatures, intrinsic networks and host factors for influenza virus A infections. NPJ Syst Biol Appl, 3, 35.

HEINZ, S., TEXARI, L., HAYES, M. G. B., URBANOWSKI, M., CHANG, M. W., GIVARKES, N., RIALDI, A., WHITE, K. M., ALBRECHT, R. A., PACHE, L., MARAZZI, I., GARCIA-SASTRE, A., SHAW, M. L. & BENNER, C. 2018. Transcription Elongation Can Affect Genome 3D Structure. Cell, 174, 1522-1536.e22.

HU, J., HU, Z., WANG, X., GU, M., GAO, Z., LIANG, Y., MA, C., LIU, X., HU, S., CHEN, S., PENG, D., JIAO, X. & LIU, X. 2018. Deep sequencing of the mouse lung transcriptome reveals distinct long non-coding RNAs expression associated with the high virulence of H5N1 avian influenza virus in mice. Virulence, 9, 1092-1111.

JOSSET, L., TCHITCHEK, N., GRALINSKI, L. E., FERRIS, M. T., EISFELD, A. J., GREEN, R. R., THOMAS, M. J., TISONCIK-GO, J., SCHROTH, G. P., KAWAOKA, Y., MANUEL DE VILLENA, F. P., BARIC, R. S., HEISE, M. T., PENG, X. & KATZE, M. G. 2014. Annotation of long non-coding RNAs expressed in collaborative cross founder mice in response to respiratory virus infection reveals a new class of interferon-stimulated transcripts. RNA Biol, 11, 875-90.

LEE, J., SHAH, M., BALLOUZ, S., CROW, M. & GILLIS, J. 2020. CoCoCoNet: conserved and comparative co-expression across a diverse set of species. Nucleic Acids Research, 48, W566-W571.

LOVE, M. I., HUBER, W. & ANDERS, S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.

SENGUPTA, S., TANG, S. Y., DEVINE, J. C., ANDERSON, S. T., NAYAK, S., ZHANG, S. L., VALENZUELA, A., FISHER, D. G., GRANT, G. R., LOPEZ, C. B. & FITZGERALD, G. A. 2019. Circadian control of lung inflammation in influenza infection. Nat Commun, 10, 4107.

STEED, A. L., CHRISTOPHI, G. P., KAIKO, G. E., SUN, L., GOODWIN, V. M., JAIN, U., ESAULOVA, E., ARTYOMOV, M. N., MORALES, D. J., HOLTZMAN, M. J., BOON, A. C. M., LENSCHOW, D. J. & STAPPENBECK, T. S. 2017. The microbial metabolite desaminotyrosine protects from influenza through type I interferon. Science, 357, 498-502.

YILDIZ, S., MAZEL-SANCHEZ, B., KANDASAMY, M., MANICASSAMY, B. & SCHMOLKE, M. 2018. Influenza A virus infection impacts systemic microbiota dynamics and causes quantitative enteric dysbiosis. Microbiome, 6, 9.

ZHAO, N., SEBASTIANO, V., MOSHKINA, N., MENA, N., HULTQUIST, J., JIMENEZ-MORALES, D., MA, Y., RIALDI, A., ALBRECHT, R., FENOUIL, R., SANCHEZ-APARICIO, M. T., AYLLON, J., RAVISANKAR, S., HADDAD, B., HO, J. S. Y., LOW, D., JIN, J., YURCHENKO, V., PRINJHA, R. K., TARAKHOVSKY, A., SQUATRITO, M., PINTO, D., ALLETTE, K., BYUN, M., SMITH, M. L., SEBRA, R., GUCCIONE, E., TUMPEY, T., KROGAN, N., GREENBAUM, B., VAN BAKEL, H., GARCIA-SASTRE, A. & MARAZZI, I. 2018. Influenza virus infection causes global RNAPII termination defects. Nat Struct Mol Biol, 25, 885-893.
