# Bacterial Genome Annotation pipeline

## Usage: `nextflow run [script] -c [*.config] -with-report [report_name].html`

## Adjustable parameters are all stored in **anot.config** file

All the configurable parameters are set in the \*.config file. Its main parameters are:

**Example of file configuration**

* Fasta file - User must provide the full path to the genome that will be annotated.
	* **params.genome** = 'GCF_000005845.2_ASM584v2_genomic.fna'
* Output folder - This parameter sets the path to store all the outputs.
	* **params.outDir** = 'out_nextflow'
* Prefix - Used to write all the output files with this prefix.
	* **params.prefix** = 'teste_1'
* Threads - Set the number of maximum threads to use.
	* **params.threads** = '4'
* Allows the use of docker containers - Must no be changed.
	* **docker.enabled** = true
* Blast Parameters - Parameters that will be used as filter.
	* **params.blast.percid** = '90'
	* **params.blast.bestHitOverhang** = '0.25'
	* **params.blast.bestHitScoreEdge** = '0.05'
	* **params.blast.eval** = '1e-10'
	* **params.blast.numDescription** = '2'
	* **params.blast.min.length** = '200'

## Further reading

A better explanation of the pipeline and its workflow is presented in the `docs/Readme.html` file. Its reading is utterly recommended in order to understand the outputs generated and how they are produced.

## Maintainer

Felipe Marques de Almeida <falmeida@aluno.unb.br>
