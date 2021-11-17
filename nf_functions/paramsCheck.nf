def paramsCheck() {
  
  /*

      Checking the prokka parameters

  */
  if (params.prokka_kingdom && !params.prokka_genetic_code) {
    println """
    ERROR!

    A minor error has occurred
      ==> User have set --prokka_kingdom but forgot --prokka_genetic_code.

    These parameters must be used together. If you change prokka defaults kingdom parameter you must set the genetic code to be used for translation.

    If in doubt with these parameters let it blank, or get more information in Prokka's documentation.

    Cheers.
    """.stripIndent()

    exit 1
  }

}
