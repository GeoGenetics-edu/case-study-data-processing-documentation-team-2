# Day 1

## Alignment with Bowtie2

First we made a symbolic link to data:
```
ln -s ~/course/data/day2/fastq/*.fq.gz .
```

Then we made a loop to remove duplicates and gzip the files:

```
for i in *.gz
do echo ${i}
j=$(basename $i .fq.gz) 
echo $j
vsearch --fastx_uniques ${i} --fastqout ${j}.vs.fq --minseqlength 30 --strand both
gzip ${j}.vs.fq 
done
```

Finally we made a loop to map the files against a database:

```
for i in *vs.fq.gz
do echo ${i}
j=$(basename $i .vs.fq.gz) 
echo $j
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U ${i} --no-unal | samtools view -bS - > ${j}.bam
done
```

# Day 2

# metaDMG

Main goal: Apply metaDMG on 5 samples efficiently.

First activate the environment and go to the working directory:

```
conda activate day1
cd ~/course/wdir/mapping
```

Then sort the five sample alignments:

```
for i in *[0-9].bam
do echo ${i}
j=$(basename $i .bam)
echo $j
samtools sort -n ${i} -@ 5 > ${j}.sort.bam
done
```

Then activate the metaDMG environment:

```
conda activate metaDMG
```

*Note: No need to run metaDMG with the ncl command. It was just to show that one could run ngsLCA directly from metaDMG.*

Run metaDMG config to create the config file to use it on all the samples. 

```
metaDMG config *.sort.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp
```

Edit the config file with `vim` and set the custom database option to `true`, the number of samples to `2` and the number of cores per sample to `2` as well. (We cannot use more than 8 cores per student at once).

Now finally run metaDMG:

```
metaDMG compute config.yaml
```

Then use metaDMG dashboard to open the visualisation:

```
metaDMG dashboard config.yaml
```

Play with the visualisation and try to see how to filter and visualise the data in a meaningful way. 

We find that terrestrial taxa generally have higher damage rates than that of aquatic taxa. Regarding Bacteria, aerobic families (e.g., Mycobacteriaceae) usually have nice damage patterns, and damage rate increases as an increase of age. However, mixed families (i.e., aerobic and anaerobic taxa), for example, Pseudomonadaceae, do not have nice damage patterns. This means that these taxa may have different taphonomic processes.

# Day 3

## Euka

First copy the folder to your own folder

```
cp -r ~/opt/vgan/share/vgan/euka_dir/ .
```

Then loop through all the fastq files and create necessary output files. Euka is basically doing an alignment to the pangenome of different clades.

```
for i in ls *.vs.fq.gz ; do j=$(basename $i .vs.fq.gz); echo $j; ~/course/data/vgan/bin/vgan euka -fq1 <(zcat ~/course/wdir/mapping/${i}) -o $j -t 5 --euka_dir ~/course/wdir/mapping/euka_dir; done
```
We didn’t go further than that I fear and we don’t have the documentation for it anymore.

## pathPhinder

We didn’t write the code for pathPhinder because it took too long to modify each piece of code and run it. 


## Population genetics

First, we activate our environment and install some R packages using mamba: 

```
conda activate day3
mamba install r-tidyverse r-ggrepel
```

Now let's set up a working folder and copy the data: 

```
mkdir popgen
cd popgen
cp ~/course/data/day3/popgen/* .
```

Calculate some basic missing data summaries with PLINK as an example of its usage:

```
plink --bfile modern_polar_mexican --missing --out modern_polar_mexican
```

The output file shows that the mexican bears have a large amount of missing data. Not sure if we can get anything out of it. 

Now we can try to run the PCA and verify this assumption using smartpca. Then we use a Rscript to plot 20 components:

```
smartpca -p modern_all.smartpca.par | tee modern_all.smartpca.log
Rscript plot_pca.R modern_all.evec label_inds.txt modern_all.pdf
```

The polar bear is really influencing strongly the first component of the PCA. The mexican bears seem to fall within the distribution of the modern black bears. They fall close to the Kenai bear and the Eastern bear. The three bears seem to be quite similar, probably because they are from the same species (or same individual)?

Now let’s rerun the analysis but without the polar bear population for better visibility.

```
smartpca -p modern_blackbear.smartpca.par | tee modern_blackbear.smartpca.log
Rscript plot_pca.R modern_blackbear.evec label_inds.txt modern_blackbear.pdf
```

Now we can visualise the repartition of the black bear species much better to try to have a view of population structure. 
The mexican black bears seem to be their own thing but they seem to be genetically closer to the Eastern population of black bears. 

Let’s venture into F-statistics to explore the data further.

We activate the R environment and install some needed packages:

```
conda activate r
mamba install bioconductor-ggtree r-ape
```

Let us use this Rscript to generate the F-statistics:

```
Rscript f_statistics.R
```

The Rscript basically goes through F4-statistics, admixture plot and pairwise Fst estimation.

F4 generates a plot with f4 on the x axis and population on the y axis. Basically F4 tests a tree with the polar bear as an outgroup, mexican bears, the test population and the eastern bears in the right. If the tree is true then the f4 statistics should be 0. If it is not 0, the scenario is more complex, like potential gene flow. It implies extra-sharing maybe between mexican and the x population. If the error bars are over the 0 line, then the tree is not rejected. However, for the case of the West and Southwest as test populations, the tree is more complicated. The value is positive, which means that the West and Southwest is genetically closer to the Eastern, and thus that the mexican bear is less closely related to the West and Southwest bears. 

