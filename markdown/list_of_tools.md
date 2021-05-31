# List of tools

These are the tools that wrapped inside bacannot. **Cite** the tools whenever you use its output.

| Analysis steps | Used software or databases |
| :------------- | :------------------------- |
| Genome assembly (if raw reads are given) |  [Flye](https://github.com/fenderglass/Flye) and [Unicycler](https://github.com/rrwick/Unicycler) |
| Identification of closest 10 NCBI Refseq genomes | [RefSeq Masher](https://github.com/phac-nml/refseq_masher) |
| Generic annotation and gene prediction | [Prokka](https://github.com/tseemann/prokka) |
| rRNA prediction | [barrnap](https://github.com/tseemann/barrnap) |
| Classification within multi-locus sequence types (STs) | [mlst](https://github.com/tseemann/mlst) |
| KEGG KO annotation and visualization | [KofamScan](https://github.com/takaram/kofam_scan) and [KEGGDecoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) |
| Methylation annotation | [Nanopolish](https://github.com/jts/nanopolish) |
| Annotation of antimicrobial (AMR) genes | [AMRFinderPlus](https://github.com/ncbi/amr/wiki), [ARGminer](https://bench.cs.vt.edu/argminer), [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/) and [RGI](https://github.com/arpcard/rgi) |
| Annotation of virulence genes |  [Victors](http://www.phidias.us/victors/) and [VFDB](http://www.mgc.ac.cn/VFs/main.htm) |
| Prophage sequences and genes annotation | [PHASTER](https://phaster.ca/) database, [Phigaro](https://github.com/bobeobibo/phigaro) and [PhySpy](https://github.com/linsalrob/PhiSpy) |
| Annotation of integrative and conjugative elements | [ICEberg](https://academic.oup.com/nar/article/47/D1/D660/5165266) |
| _In silico_ detection of plasmids | [Plasmidfinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/) and [Platon](https://github.com/oschwengers/platon) |
| Prediction and visualization of genomic islands | [IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath) and [gff-toolbox](https://github.com/fmalmeida/gff-toolbox) |
| Focused detection of insertion sequences | [digIS](https://github.com/janka2012/digIS) |
| Merge of annotation results | [bedtools](https://bedtools.readthedocs.io/en/latest/) |
| Renderization of results in a Genome Browser | [JBrowse](http://jbrowse.org/) |
| Renderization of automatic reports and shiny app for results interrogation | [R Markdown](https://rmarkdown.rstudio.com/) and [Shiny](https://shiny.rstudio.com/) |
