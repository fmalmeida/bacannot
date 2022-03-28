# Output files

Here, using the results produced in the [quickstart section](quickstart.md#), we give users a glimpse over the main outputs produced by bacannot. The command used in the quickstart wrote the results under the `_ANNOTATION` directory.

!!! note

    Please take note that the pipeline uses the directory set with the `--output` parameter as a storage place in which it will create a folder for each sample using its `id`. Therefore the the same `--output` can be used for different annotations.

## Directory tree

After a successful execution, you will have something like this:

```bash

# Directory tree from the running dir
.
├── _ANNOTATION
|   └── ecoli_ref.fna
│   └── ecoli
│       ├── assembly          # Assembly files (when raw reads are given)
│       ├── annotation        # Prokka annotation files
│       ├── antiSMASH         # antiSMASH secondary annotation files
│       ├── digIS             # Insertion sequences predicted with digIS
|       ├── gbk               # Gbk files produced from the resulting GFF
|       ├── gffs              # A copy of the main GFF files produced during the annotation
|       ├── genomic_islands   # Genomic Islands predicted with IslandPath-DIMOB
|       ├── ICEs              # Results from ICEberg database annotation
|       ├── jbrowse           # The files that set up the JBrowse genome browser
|       ├── KOfamscan         # Results from annotation with KEGG database
|       ├── methylations      # Methylated sites predicted with Nanopolish (if fast5 is given)
|       ├── MLST              # MLST results with mlst pipeline
|       ├── plasmids          # Plasmid annotation results from Platon and Plasmidfinder
|       ├── prophages         # Prophage annotation results from PhiSpy, Phigaro and PHAST
|       ├── refseq_masher     # Closest NCBI Resfseq genomes identified with refseq_masher
|       ├── report_files      # Annotation reports in HTML format
|       ├── resistance        # AMR annotation results from ARGminer, AMRFinderPlus, RGI and Resfinder
|       ├── rRNA              # barrnap annotation results
|       ├── SequenceServerDBs # SequenceServer pre-formatted databases to be used with SequenceServer blast application
|       ├── SQLdb             # The SQLdb of the annotation used by the shiny server for rapid parsing
|       ├── tools_versioning  # Versions of tools and databases used (whenever available)
|       ├── virulence         # Virulence genes annotation results from Victors and VFDB databases
|       └── run_server.sh     # The shiny parser runner that enables a rapid and simple exploration of the results (see below)
```

## KEGG KO annotation heatmap

Using both [KofamScan](https://github.com/takaram/kofam_scan) and [KEGGDecoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder), bacannot is capable of annotating KOs and plotting a heatmap of detected pathways as exemplified below.

!!! tip ""

	Click on the image to zoom it! :)

<center>
  <img src="../images/ecoli_kegg-decoder_heatmap-static.svg" width="100%">
</center>

## Bacannot automatic reports

Bacannot will use [R Markdown](https://rmarkdown.rstudio.com/) to produce automatic annotation reports. To date, the available reports are:

* [Report of general annotation features](https://fmalmeida.github.io/reports/report_general.html)
* [Report of Antimicrobial resistance (AMR) genes annotation](https://fmalmeida.github.io/reports/report_resistance.html)
* [Report of virulence genes annotation](https://fmalmeida.github.io/reports/report_virulence.html)
* [Report of mobile genetic elements annotation](https://fmalmeida.github.io/reports/report_MGEs.html)
    * Including plasmids, prophages, ICEs and genomic islands.
* Report of user's custom db annotations.
    * The quickstart does not produce an example, however, the report is similar to the ICEberg section in the MGE example report.
    * See [custom-db reference page](custom-db.md#)
* [Report of antiSMASH annotation](https://docs.antismash.secondarymetabolites.org/understanding_output/)
    * The annotation report is provided by the antiSMASH tool

## Genome Browser

With aid of [JBrowse](http://jbrowse.org/), Bacannot already give users a totally customised and redered Genome Browser for exploration of annotation results.

<center>
  <img src="../images/jbrowse.png" width="85%">
</center>

!!! warning

    The JBrowse wrapper in the shiny server is not capable of displaying the GC content and methylation plots when available. It can only display the simpler tracks. If the user wants to visualise and interrogate the GC or methylation tracks it must open the JBrowse outside from the shiny server. For that, two options are available:
      
    * You can navigate to the `jbrowse` directory under your sample's output folder and simply execute `http-server`. This command can be found at: https://www.npmjs.com/package/http-server
    * Or, you can download the `JBrowse Desktop app <https://jbrowse.org/docs/jbrowse_desktop.html) and, from inside the app, select the folder `jbrowse/data` that is available in your sample's output directory.

In order to provide an integrative solution, the genome browser is already packed inside the shiny app that can be launched with the `run_server.sh` script that loads the server docker image (See below at Bacannot shiny parser).

## Bacannot shiny parser

<center>
  <img src="../images/bacannot_shiny.gif" width="70%">
</center>

The bacannot shiny server is basically a wrapper of the main outputs of the pipeline that is packed up in a docker image called `fmalmeida/bacannot:server`. This server is triggered by going under the results folder, in our quickstart case for instance, the `_ANNOTATION/ecoli` folder, and executing the command:

```bash
# Trigger the server
./run_server.sh -s

# This will open the pipeline in localhost:3838
# log message:
The server has started in: http://localhost:3838/
When finished, run the command:
	    docker rm -f ServerBacannot

# To stop the server you just need to execute
docker rm -f ServerBacannot
```

###  Server homepage

In the first page of the shiny app, the main HTML reports and the **JBrowse genome browser** are indexed as url links for quick opening (See the image below).

<center>
  <img src="../images/bacannot_server_home.png" width="85%">
</center>

### Server SQLdb parser

In the second page, the SQL database (SQLdb) produced in the pipeline is used to provide a rapid and simple way to query and filter the genome annotation.

!!! note

    The SQLdb parser contains a set of features that enables users to filter the annotation following their desires. It is possible to filter based on `contigs`, `sources`, `start`, `end`, `strand` and more.

    Additionally, it accepts as input a file of patterns. These patterns are used to filter the annotation based on the values available in the attributes column of the GFF (9th column).
    
    Any value available in this column can be used as filters, the only requirement is to write each pattern in one line, exactly as it is found in the annotation result. For example, it can be used to select only a few genes based on their IDs.

<center>
  <img src="../images/bacannot_server_sqldb.png" width="85%">
</center>

### Server BLAST (for intersection) app

In the its third page, the server provides a simple way to BLAST the genome with new queries and to automatically identify intersections between the blast results and the the main annotation.

<center>
  <img src="../images/bacannot_server_blast.png" width="85%">
</center>

### Server BLAST (SequenceServer) app

In its the last page, the server provides an implementation of [SequenceServer](https://sequenceserver.com/) which allows users to BLAST their samples and visualise the alignments produced.

<center>
  <img src="../images/bacannot_server_blast_sequenceserver.png" width="85%">
</center>
