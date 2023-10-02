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
