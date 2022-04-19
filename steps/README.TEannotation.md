 # Cactoblastis genome sequencing project
 Genome sequencing project of species of the Cactoblastis genus

 ## Index
 + 00.RawData
 + 01.QC
 + 02.Assembly
 + 03.RNA-mapping
 + 04.TE-annotation
 + 05.Annotation
 + 06.Citations

 ## 00.RawData

 ## 01.QC

 ## 02.Assembly

 ## 03.RNA-mapping

 ## 04.TE-annotation
 Transposable Elements (TEs) annotation was performed following the pipeline described in [Bell et al. (2021)](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13489). We applied this approach with four genomes of cactophilic *Drosophila* species of the *buzzatii* cluster (*repleta* group) that we recently sequenced and reported in [Moreyra et al. (2022)]().
 We describe here each stage of TE annotation, most of which were adapted from [FasTE](https://github.com/ellenbell/FasTE).

 The script employed in each stage was ran with Linux Ubuntu (v18.04.5 LTS Bionic Beaver), 36 cores, 512GB RAM on each genome (size ranging 166-190MB). 

 ### 04.00. EDTA: TE library generation
 We used the Extensive de-novo TE Annotator (EDTA) v1.9.4 [(Ou et al. 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y) to generate a whole-genome de-novo TE annotation.

 ~~~~
 #!/bin/bash
 
 WD=/home/nmoreyra/Data/07_TE-annotation/EDTA
 LIB=nocuratedlib
 LIB="${2:-$LIBDEFAULT}"

 for species in Dato Dbrb DkoeA DkoeB; do
    cd $WD/${species}_${LIB}
    echo "#####${species}"
    EDTA.pl --genome ${species}_genome.fasta --cds ${species}_CDS.fasta --exclude ${species}_annotations.norepeats.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 36 &> log.${species}.out
 done
 ~~~~

 On this system with this genome, EDTA ran about XX hours for the four genomes (~160-190 Mb).

 ### 04.01 DeepTE: TE Classification
 We employed DeepTE [(Yan et al. 2020)](https://academic.oup.com/bioinformatics/article/36/15/4269/5838183?login=true) to classify transposons with unknown classification.

 ~~~~
 #!/bin/bash

 species=$1
 LIBDEFAULT=nocuratedlib
 LIB="${2:-$LIBDEFAULT}"
 WD=/home/nmoreyra/Data/07_TE-annotation/DeepTE
 EDTA=/home/nmoreyra/Data/07_TE-annotation/EDTA
 MDIR=/home/nmoreyra/Software/DeepTE/Model_dir/Metazoans_model

 cd $WD/${species}
 mkdir -p $WD/${species}/${LIB} $WD/${species}/${LIB}/tmp
 echo "#####${species}"
 DeepTE.py -d $WD/${species}/${LIB}/tmp -o $WD/${species}/${LIB} -i ${EDTA}/${species}_${LIB}/${species}_genome.fasta.mod.EDTA. TElib.fa -sp M -m_dir ${MDIR}
 
 #Additional step to clean up FASTA headers obtained after EDTA and DeepTE
 sed -e 's/\(#\).*\(__\)/\1\2/'  opt_DeepTE.fasta > opt_DeepTE.headercleaned.fasta

 echo ">>>${species}_$LIB has finished"
 ~~~~

 DeepTE ran in about XX hours to classify four TE libraries (~1.8-2.2 Mb).
 
 ### 04.02 RepeatMasker
 Now we use [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) (RM) v4.1.1 to screen for TEs using the generated de novo TE library. 

 ~~~~
 #!/bin/bash
 WD=/home/nmoreyra/Data/07_TE-annotation/RM
 EDTA=/home/nmoreyra/Data/07_TE-annotation/EDTA #/Dbrb_nocuratedlib/Dbrb_genome.fasta
 LIBDEFAULT=nocuratedlib
 LIB="${2:-$LIBDEFAULT}"
 
 for species in Dato Dbrb DkoeA DkoeB; do
   cd $WD/${species}/${LIB}
   echo "#####${species}"
   RepeatMasker $EDTA/${species}_${LIB}/${species}_genome.fasta -e rmblast -pa 36 -s -lib $WD/../DeepTE/${species$
   
   ##remove lines with an asterisk in them (repeats that overlap with one or more other hits that have a higher score)
   awk '!/\*/' ${species}_genome.fasta.out > ${species}_genome.fasta.noasterisk.out
   
   ##Remove -int notation from the TE name
   sed 's/-int//' ${species}_genome.fasta.noasterisk.out > ${species}_genome.fasta.noasterisk.tidy.out
 done
~~~~

 ### 04.03 RM_TRIPS
 Parsing RepeatMasker Output with [RM_TRIPS](https://github.com/clbutler/RM_TRIPS).

 ~~~~
 #!/bin/bash
 
 WD=/home/nmoreyra/Data/07_TE-annotation/TRIPS
 RM=/home/nmoreyra/Data/07_TE-annotation/RM
 LIBDEFAULT=nocuratedlib
 LIB="${2:-$LIBDEFAULT}"

 for species in Dato Dbrb DkoeA DkoeB; do
   cd $WD/${species}/${LIB}
   echo "#####${species}"
   #nohup time RepeatMasker $EDTA/${species}_${LIB}/${species}_genome.fasta -e rmblast -pa 36 -s -lib $WD/../DeepTE/${species}/${LIB}/opt_DeepTE.headercleaned.fasta -dir . &> RM.${species}.${LIB}.log
   ### set up inputs
   i=/home/nmoreyra/Data/07_TE-annotation/RM/${species}/$LIB #directory where .out file is located
   j=${species}_genome.fasta.noasterisk.tidy.out #set name of file
   k=/home/nmoreyra/Data/07_TE-annotation/DeepTE/${species}/$LIB #directory where the repeatmasker library is found (.lib/fasta file)
   l=opt_DeepTE.headercleaned.fasta #set name of .lib file
   #run RM_TRIPS
   nohup time Rscripts RM_TRIPS.args.R $i $j $k $l &> ${species}.${LIB}.log
 done
 ~~~~
 ## 05.Annotation

 ## 06.Citations
 + Bell, Ellen A., Christopher L. Butler, Claudio Oliveira, Sarah Marburger, Levi Yant, and Martin I. Taylor. "Transposable element annotation in non‐model species: The benefits of species‐specific repeat libraries using semi‐automated EDTA and DeepTE de novo pipelines." Molecular Ecology Resources 22, no. 2 (2022): 823-833.
 + Moreyra et al. (2022)
 + Ou, Shujun, Weija Su, Yi Liao, Kapeel Chougule, Jireh RA Agda, Adam J. Hellinga, Carlos Santiago Blanco Lugo et al. "Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline." Genome biology 20, no. 1 (2019): 1-18.
 + Yan, Haidong, Aureliano Bombarely, and Song Li. "DeepTE: a computational method for de novo classification of transposons with convolutional neural network." Bioinformatics 36, no. 15 (2020): 4269-4275.

 

