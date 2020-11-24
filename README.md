# Práctico SNPs

Introducción

En el siguiente práctico ejecutaremos el pipeline de [GATK](https://gatk.broadinstitute.org/hc/en-us).
Luego anotaremos las variantes y buscaremos el resultado en bases de datos internacionales.

Debe crear una carpeta en su cuenta:

      mkdir GATK

Entre a su carpeta y descargue GATK.


Deberan descargar [GATK 4.1.9](https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip):


      wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
     
     
Realice un unzip de su carpeta:

      unzip gatk-4.1.9.0.zip
     
Quedará una carpeta gatk-4.1.9.0, revise su carpeta y ejecutemos un test, 
pero antes necesitaremos cargar java 8 en nuestro ambiente para probar el programa:

Revisemos las versiones de Java:
  
     module avail
     
 
Existen varias versiones de Java, cargaremos Java:

    ml Java/1.8










