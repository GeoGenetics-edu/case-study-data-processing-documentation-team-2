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


[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)


#Team 2, day 2, Tuesday

Main goal: Change script to deal with 5 samples instead of 1

First activate environment and go to working directory

```
conda activate day1
cd ~/course/wdir/mapping
```

Then sort the five sample alignments

```
for i in *[0-9].bam
do echo ${i}
j=$(basename $i .bam)
echo $j
samtools sort -n ${i} -@ 5 > ${j}.sort.bam
done
```

Then activate the metaDMG environment

```
conda activate metaDMG
```

Note: No need to run the metaDMG command with ncl

Run metaDMG config to create the config file to use it on all the samples. 

```
metaDMG config *.sort.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp
```

Edit the config file with vim and set the custom database option to true, the number of samples to 2 and the number of cores per sample to 2 as well. (We cannot use more than 8 cores per student at once).

Now finally run metaDMG

```
metaDMG compute config.yaml
```

Then use metaDMG dashboard to open the visualisation.

```
metaDMG dashboard config.yaml
```

Play with the visualisation and try to see how to filter and visualise the data in a meaningful way. 
