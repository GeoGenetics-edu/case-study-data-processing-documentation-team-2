[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)
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

## Data processing and visualisation of the metaDMG results

1. Convert the output of metaDMG to a merged csv file.

```
cd ~/course/wdir/mapping
metaDMG convert --output metaDMGresults.csv --add-fit-predictions
```

2. Create a new working directory

```
mkdir ~/course/wdir/mapping/plots
```

3. Activate r environment and R.

```
conda activate r
R
```
We actually ran all R scripts in R Studio.

4. Activate libraries.

```
library(tidyverse) 
library(reshape2)
library(vegan)
library(rioja)
library(ggplot2)
library(dplyr)
library(gghighlight)
```

5. Set working directory.

```
setwd("~/course/wdir/mapping/plots/")
```

6. Import data and metadata

```
df <- read_csv("~/course/wdir/mapping/metaDMGresults.csv")
metaDATA <- read.delim("~/course/data/shared/metadata/metadata.tsv")
```

Then we checked the count number of columns in the dataframe 'df', which is 187.

7. Replace header "sample_name" with "sample"

```
colnames(metaDATA)[colnames(metaDATA) == "sample_name"] <- "sample"
colnames(metaDATA)[colnames(metaDATA) == "years_bp"] <- "YearsBP"
```

8. Merge the metadata and dataframe by "sample"

```
dt <- merge(df, metaDATA, by = "sample")
```

9. Then we checked if new columns have been added to the dt dataframe,

```
ncol(df) < ncol(dt)
```

and the result is [1] TRUE.

10. Set the parameters that to filter the data

```
DamMin2 = 0.00
MapSig2 = 0
MinRead2 = 100
MinLength = 35
```

The reason we used DamMin2=0 and MapSig2=0 is that we would like to keep all potential damage information at very beginning. We don't consider the results with read count less than 100. Also the fragment length should longer than 35 bp.

11. Subset the table by Viridiplantae

```
dt2 <- dt %>% filter(MAP_damage > DamMin2, N_reads >= MinRead2, mean_L > MinLength, MAP_significance  > MapSig2,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))
```

12. Plot the results of Viridiplantae

```
pdf(file = "aeCourse.DNAdamageModelJitterPlot.pdf", width = 8, height = 4)
ggplot() +
  geom_jitter(data = dt2, aes(x=as.numeric(YearsBP), y=MAP_damage, size = N_reads), alpha =0.5) +
  gghighlight(N_reads > 500) +
  xlab("Years BP") +
  ylab("DNA damage") +
  labs(color = "Values for taxa with \n>500 reads", size = "Number of reads")
dev.off()
```

![Figure1](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team-2/assets/146089734/109711e8-e808-4ad4-b9a9-301aae1f4550)

We can see damage rate generally increases with age, which makes sense.

13. Plot plant taxa, highlight taxa with more than 500 reads and add the min, max and median.

```
pdf(file = "aeCourse.DNAdamageLRJitterPlot.pdf", width = 8, height = 4)
ggplot() +
  geom_jitter(data = dt2, aes(x=as.numeric(YearsBP), y=MAP_damage, size = MAP_significance), alpha =0.5) +
  gghighlight(N_reads > 500) +
  xlab("Years BP")+
  ylab("DNA damage") +
  labs(color = "Values for taxa with \n>500 reads", size = "Significance \nfor Taxa with >500 reads")
dev.off()
```

![Figure2](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team-2/assets/146089734/bc9b4ff8-60be-4831-8061-2e284bdb2d4d)

Significance shows the same pattern. In addition, the damage in the older samples seems to be more convincing (significant) than that in the younger samples. 

14. Create filtered table for DNA damage model

```
filtered_data <- dt2 %>% filter(N_reads >= 500)
```

15. Set the parameters that to further filter the data

```
MapSig3 = 3
MinRead3 = 100
MinLength3 = 35
```

This time we increase the significance to 3.

16. Further subset the table

```
filtered_data_viridiplantae <- filtered_data %>% filter(N_reads >= MinReads3, mean_L > MinLength3, MAP_significance > MapSig3,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", YearsBP))
```
Then we counted the number of unique plant taxa, which is 78.

17. Prepare a table for downstream plot and data wrangling of the plants.

```
data_wide_plants <- dcast(filtered_data_viridiplantae, tax_name ~ YearsBP, value.var="N_reads", fun.aggregate = sum)
n <- ncol(data_wide_plants)
b2 <- data_wide_plants[,2:n]
rownames(b2) <- data_wide_plants$tax_name
b2[is.na(b2)] <- 0 #set all NAs as zeros
```

18. Calculate the precentage (relative abundance) of each taxon in each sample.

```
i=ncol(b2)
b3=as.matrix(b2[,seq(1,i)])  
b4 <- prop.table(data.matrix(b3), margin=2)*100
colSums(prop.table(b4, margin=2)*100)
```

19. Plot the strat.plot

```
b5 <- t(b4)
z <- as.numeric(rownames(b5)) # depth/depth
pdf(file = "aeCourse.Stratplot_Plants_area.pdf", width = 15, height = 5)
pdf(file = "aeCourse.Stratplot_Plants_area1.pdf", width = 40, height = 10)
strat.plot(b5, y.rev=TRUE, plot.line=TRUE, plot.poly=TRUE,
           plot.bar=FALSE, lwd.bar=10, sep.bar=F, scale.percent=T,
           xSpace=0.0001, x.pc.lab=TRUE, x.pc.omit0=TRUE, srt.xlabel=90,
           las=2, exag=TRUE, exag.mult=5, ylabel = "years BP",  yvar = z)
dev.off()
```

![Figure3](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team-2/assets/146089734/b5d6c43d-e678-48be-90ba-7f94c1ef8d17)
Hum... Not very nice. Let's plot the heatmap.

20. Plot the heatmap

```
y <- ncol(b5)
b6 <- melt(b5[,1:y])
sapply(b6, class)
colnames(b6) <- c("YearsBP","Taxa", "percentage")
p1 <- ggplot(b6, aes(y=Taxa, x=YearsBP, fill=percentage)) + 
  geom_tile(colour="lightgrey") +
  theme_minimal() + scale_fill_gradient(low="white", high="darkgreen")+
scale_y_discrete(limits=rev)
p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + 
  ggtitle("Percentage of taxa plotted as heatmap") +
  xlab("YearsBP") + ylab("Taxa name") + labs(fill = "percentage %")
```

![Percentage of taxa plotted as heatmap](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team-2/assets/146089734/bf9ba92f-d3f1-4694-a1cc-3edb43bffda0)

Much better! We can see that plant community in this cave changed fast in the past 6000 years. The two younger samples have higher richness than the older samples. In addition, plant composition seems strange, including many tree and shrub taxa (e.g. Alnus, Betula, Populus, Salix), even aquatic plant taxa (e.g. Hippuris). Therefore, the sampling site was probably close to the cave entrance?

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

