process mongoDB {
  validExitStatus 0,48
  publishDir "${params.outdir}/${prefix}/mongoDB", mode: 'copy'
  label 'main'

  input:
  tuple val(prefix), file(gff)

  output:
  tuple val(prefix), file("data") // Get all files created in the MongoDB server
  tuple val(prefix), file("${prefix}.json") // Save JSON file for future uses

  script:
  """
  # Activate Python Env
  source activate mongoDB

  # Create required directories
  mkdir data data/db

  # Create mongoDB server
  mongod --dbpath data/db --logpath data/mongo_log.log --fork

  # Convert GFF to JSON
  gff2json.py -i $gff -o ${prefix}.json

  # Create MongoDB Collections
  mongoDB_parse_JSON.py -i ${prefix}.json -n ${prefix}

  # Close server
  mongod --dbpath data/db --logpath --shutdown 2> /tmp/mongo.tmp
  """
}
