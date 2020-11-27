# Práctico SNPs

Introducción

En el siguiente práctico ejecutaremos el pipeline de [GATK](https://gatk.broadinstitute.org/hc/en-us).
Luego anotaremos las variantes y buscaremos el resultado en bases de datos internacionales.

#### Objetivos

1. Demostrar la capacidad de instalar un programa en HPC desde su carpeta local sin permiso de root.
2. Poder ejecutar el pipeline y la generación de scripts `SBATCH` por si solos.
3. Montar la capacidad de buscar variantes en HPC.
4. Buscar Variantes utilizando GATK.
5. Anotar las variantes con ANNOVAR.



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

Podemos monitorear el estado de nuestro job con el comando watch y squeue:

    watch squeue

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
Para nuestro beneficio estos requermientos vienen en un archivo de configuracion yaml.

Nos vamos a:

    cd GATK
    cd gatk-4.1.9.0
    
Levantamos el modulo de java:

    ml Java/1.8
    
 y crearemos un ambiente gatk:
 
    conda env create python=3.6.2 -f gatkcondaenv.yml

Con esto hemos creado un nuevo ambiente python 3.6.2, puede revisar sus ambientes en `~/anaconda3/envs/` o con ```conda environments```.

Para activar el ambiente ejecutaremos

    conda activate gatk
      
el prompt debería cambiar a:
      
    (gatk)user@server:~$

podemos ver que nuestro python ahora cambio de dirección de nuevo:

    which python
      
ahora volvemos a la carpeta de GATK:

    cd GATK

## GATK Pipeline Magic

##### Basado en las recomendaciones de [Ricardo Palma](http://www.cmm.uchile.cl/?cmm_people=ricardo-palma)


### Paso 1. Filtrado:

Realizaremos el filtrado con [BBDUK](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) “””aka Magical Mathematics of Console Gymnastics ””” 

    bbduk.sh -Xmx3g in=R1.fq in2=R2.fq ref=adaptor_file.fa mm=f \\
    rcomp=f out=clean_R1 out2=clean_R2 threads=20 minlen=read_min_len \\
    qtrim=lr trimq=20 ktrim=r k=21 mink=9 hdist=1 tpe tbo overwrite=true

Parámetros:

read_min_len: Minimal size of read length before trimming. (depends on sample 
profiles)

### Paso 2, Mapping [BWA](http://bio-bwa.sourceforge.net/)

    bwa index reference.fasta
    bwa mem -t 20 -o output.sam clean_R1.fq clean_R2.fq
    samtools view -bS output.sam > output.bam


### Paso 3, [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037224932-MarkDuplicatesSpark)

Marcamos los duplicados  y ordenamos el BAM

    gatk MarkDuplicatesSpark -I output.bam -O marked_duplicates_sorted.bam

#### Opcional

Si es que no llegara a funcionar podemos utilizar siempre el método antiguo

    java -jar picard.jar MarkDuplicates I= output.bam O=marked_duplicates.bam M=marked_dup_metrics.txt

    samtools sort -T temp_sorted.bam -o marked_duplicates_sorted.bam marked_dup_metrics.txt


### Paso 4 [Haplotypecalller](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller): 

Variant Calling Per-Sample


    gatk –java-options “-Xmx12G” HaplotypeCaller -R reference.fasta -I marked_duplicates_sorted.bam -O output.g.vcf.gz -ERC GVCF

### Paso 5 

VariantRecalibrator, and [ApplyVQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants)
-either-with-VQSR-or-by-hard-filtering 

#### First, split the g_vcf. Files into SNPs and INDELs:

    gatk SelectVariants -V output.g.vcf.gz -select-type SNP -O only_snps.vcf.gz
    gatk SelectVariants -V output.g.vcf.gz -select-type INDEL -O only_indels.vcf.gz

Now apply the “Hard-filtering” for Variants.


##### Hard-Filtering for SNPs

    gatk VariantFiltration -V only_snps.vcf.gz -filter "QD < 2.0" --filter-name "QD2" \\
    -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" \\
    -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" \\
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" \\
    --filter-name "ReadPosRankSum-8" -O snps_filtered.vcf.gz


##### Hard-Filtering for INDELs

    gatk VariantFiltration V only_indels.vcf.gz -filter "QD < 2.0" --filter-name "QD2" \\
    -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" \\
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O indels_filtered.vcf.gz


    If a variant Pass all filters it will be tag as PASS in the “INFO” field 
    of the VCF, if not the name of the fault filter will be displayed.

### Last Step

Extraer todas las variantes que pasaron el filtro `PASS` en GATK.


###### _Congrats, you skipped 6 months navigating GATK forums. Ricardo Palma_ 



      






      

       
      








