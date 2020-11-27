# Práctico SNPs

#### Introducción

En el siguiente práctico ejecutaremos el pipeline de [GATK](https://gatk.broadinstitute.org/hc/en-us).
Luego anotaremos las variantes y buscaremos el resultado en bases de datos internacionales.

#### Objetivos

1. Demostrar la capacidad de instalar un programa en HPC desde su carpeta local sin permiso de root.
2. Poder ejecutar el pipeline y la generación de scripts `SBATCH` por si solos.
3. Montar la capacidad de buscar variantes en HPC.
4. Buscar Variantes utilizando GATK ([artículo](https://genome.cshlp.org/content/20/9/1297.short)).
5. Anotar las variantes con [ANNOVAR](http://wannovar.wglab.org/) ([artículo](https://www.nature.com/articles/nprot.2015.105)).
6. Interpretar los resultados.

#### Fechas

1. La entrega será el 14 de diciembre del 2020 (no podemos pasarnos de esa fecha), sin embargo,
2. Las cuentas estarán disponibles hasta el 16 de enero del 2021, por si les interesa utilizar el clúster.



## Preliminares

El ejercicio comienza desde cero, es decir, generaremos un ambiente en el [NLHPC](www.nlhpc.cl) por nuestra cuenta para ejecutar GATK.
Para esto, deberemos preparar todo un ambiente. 

Conéctese al NLHPC utilizando su usuario y contraseña y en su `HOME`... 

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
     
 
 Podemos ver que existen varias versiones de Java, cargaremos la siguiente:

    ml Java/1.8


Crearemos un ambiente python para instalar algunas dependencias necesarias para ejecutar GATK:

    cd ..
    mkdir python
    cd python


Para esto descargaremos miniconda, iremos a [conda](https://docs.conda.io/en/latest/miniconda.html) y descargaremos la última versión en este caso 3.8:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh



Crearemos un job para instalar anaconda:

    vim install.sbatch
      
 escribimos lo siguiente (deberá reemplazar su usuario y mail en el campo mail, para revisar que le llega la secuencia):
 ```bash
 #!/bin/sh

 #SBATCH --job-name=miniconda
 #SBATCH --output=mini.out
 #SBATCH --error=mini.err
 #SBATCH --mail-user=user@mail.com
 #SBATCH --nodes=1
 #SBATCH --mail-type=ALL
 sh Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/anaconda3
```

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
El primer caveat es que la versión requerida para este GATK, es python 3.6.2 descrita en el [README.md](https://github.com/broadinstitute/gatk) de GATK.
Por lo que vamos a generar un ambiente virtual de nuestro miniconda que contenga solo los requerimientos de GATK.
Para nuestro beneficio estos requermientos vienen en un archivo de configuracion yaml.

Nos vamos a:

    cd GATK
    cd gatk-4.1.9.0
    
Levantamos el módulo de java:

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
    
Todavía nos queda instalar algunas herramientas bioinformáticas que utilizaremos en el práctico.

    mkdir bin
    cd bin

Instalaremos [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

    wget https://sourceforge.net/projects/bbmap/files/latest/download

Podemos ver que al descargar nos quedará un archivo `download`.

si ejecutamos un file:
    
    file download
    
veremos que es `download: gzip compressed data, was "BBMap_38.87.tar"`, es decir un archivo `.tar.gz`

    tar xvfz download
    rm download
    cd bbmap

revisaremos si funciona bbduk:
    
    ./bbduk.sh

Ahora volvemos a nuestro `$HOME` y crearemos la carpeta `READS`:

    cd
    mkdir reads

Dentro de la carpeta reads descargaremos los siguientes archivos del proyecto 1000 genomas:

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

Luego extraeremos los reads:

     gzip -d SRR062634_1.filt.fastq.gz
     gzip -d SRR062634_2.filt.fastq.gz

Haremos un enlace simbólico a los reads:

    ln -s SRR062634_1.filt.fastq.gz R1.fq
    ln -s SRR062634_2.filt.fastq.gz R2.fq


## GATK Pipeline Magic

##### Basado en las recomendaciones de [Ricardo Palma](http://www.cmm.uchile.cl/?cmm_people=ricardo-palma)


### Paso 1. Filtrado:

Para cada comando o conjunto de comandos deberá realizar un script `SBATCH` para ejecutar sus tareas.
Realizaré solo el primer script como instructivo:

Realizaremos el filtrado con [BBDUK](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) “””aka Magical Mathematics of Console Gymnastics ””” 

    bbduk.sh -Xmx3g in=R1.fq in2=R2.fq ref=adaptor_file.fa mm=f \\
    rcomp=f out=clean_R1 out2=clean_R2 threads=20 minlen=read_min_len \\
    qtrim=lr trimq=20 ktrim=r k=21 mink=9 hdist=1 tpe tbo overwrite=true

Parámetros:

read_min_len: Minimal size of read length before trimming. (depends on sample 
profiles)

Script `SBATCH`:

Crearemos un archivo 
    
    vim runbb.sbatch

y escribimos el siguiente contenido:

```sh
#!/bin/sh

#SBATCH --job-name=bbduk ### ASIGNAMOS EL NOMBRE: bbduk
#SBATCH --output=bbduk.out ### asignamos un archivo que nos guarde la salida o el STDOUT del programa
#SBATCH --error=mini.err ### asignamos un archivo que nos guarde el error o el STDERR del programa
#SBATCH --mail-user=user@mail.com ### asignamos un mail para que nos envié el status del job
#SBATCH --nodes=1   ### forzamos a que esto se ejecute en un solo nodo
#SBATCH --mail-type=ALL ### queremos que nos lleguen todos los status
#SBATCH --mem=4G ### en el flag -Xmx3g de bbduk estamos asignando 3 gigabytes de ram por lo que asignamos 4G de ram ante un inesperado peak.
#SBATCH --cpus-per-task=20 ### en el flag threads=20 de bbduk estamos asignando 20 cpus por lo que tenemos que configurarlo en el scheduler.
 
### INDICAMOS AL NODO DONDE ESTA BBDUK
export PATH=~/bin/bbmap:$PATH

#### Ejecutamos bbduk

bbduk.sh -Xmx3g in=R1.fq in2=R2.fq ref=adaptors.fa mm=f rcomp=f out=clean_R1 out2=clean_R2 threads=20 minlen=read_min_len qtrim=lr trimq=20 ktrim=r k=21 mink=9 hdist=1 tpe tbo overwrite=true

```

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



      






      

       
      








