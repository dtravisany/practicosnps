# Práctico SNPs

Introducción

En el siguiente práctico ejecutaremos el pipeline de [GATK](https://gatk.broadinstitute.org/hc/en-us).
Luego anotaremos las variantes y buscaremos el resultado en bases de datos internacionales.


## Preliminares

El ejercicio comienza desde cero, es decir, generaremos un ambiente en el NLHPC por nuestra cuenta para ejecutar GATK.
Para esto, deberemos preparar todo un ambiente.

Comenzaremos creando una carpeta inicial:

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


Crearemos un ambiente python para instalar algunas dependencias necesarias para ejecutar GATK:

    cd ..
    mkdir python
    cd python


Para esto descargaremos miniconda, iremos a [conda](https://docs.conda.io/en/latest/miniconda.html) y descargaremos la última versión en este caso 3.8:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh



Crearemos un job para instalar anaconda:

    vim install.sbatch
      
 escribimos lo siguiente:
 
    #!/bin/sh

    #SBATCH --job-name=miniconda
    #SBATCH --output=mini.out
    #SBATCH --error=mini.err
    #SBATCH --mail-user=user@mail.com
    #SBATCH --nodes=1
    #SBATCH --mail-type=ALL
    sh Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/anaconda3

Guardamos y ejecutamos nuestro job:

    sbatch install.sbatch
      
y esperamos la instalación, luego de que se haya completado la instalación (aprox 10 min).
Si lo ejecutó por defecto, su instalación debería quedar en el path: `/home/courses/student[07-10]/anaconda3`.

Cuando la instalación se haya completado, inicializaremos conda para que considere 
nuestro "brand new" python como el python por defecto:

      conda init

Nos saldremos de NLHPC y volveremos a hacer login

      exit
      ssh student[07-10]@server

Ahora deberiamos ver en nuestro prompt algo como esto:

      (base)user@server:~$
      
Ya tenemos nuestro python local corriendo y listo para instalarle paquetes.
El primer caveat es que la versión requerida para este GATK, es python 3.6.2 descrita en el README.md de GATK.
Por lo que vamos a generar un ambiente virtual de nuestro miniconda que contenga solo los requerimientos de GATK.
De esta manera no interferimos con los paquetes ya instalados, para esto haremos un 
      
      conda create -n GATK python=3.6.2
      
Para activar el ambiente ejecutaremos

      conda activate GATK
      
el prompt debería cambiar a:
      
      (GATK)user@server:~$

podemos ver que nuestro python ahora cambio de dirección de nuevo:

      which python
      


Con esto hemos creado un nuevo ambiente python 3.6.2, puede revisar sus ambientes en:

      cd ~/anaconda3/envs/
      




Para esto utilizaremos el paquete `virtualenv` que deberemos descargar mediante pip:

      pip install virtualenv



      

       
      








