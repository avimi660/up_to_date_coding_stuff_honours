```
## up_to_date_coding_stuff_honours


# this is the initial completed bash script that we worked on to pull gene names from the possum genome based off the human genome annotations 
## take note that I am using specific versions of the human reference genome and the possum genome 

no_mito_genes=`wc -l mito_genes.txt | awk '{ print $1 }'`

for mito_line in `seq 2 1 $no_mito_genes`;
      do gene_search_term=`head -n $mito_line mito_genes.txt | tail -n 1 | cut -f 4`;
      grep "ID=gene-"$gene_search_term";" /nesi/nobackup/uoo03398/michael/possum_genome_master/GCF_011100635.1_mTriVul1.pri_genomic.gff | grep $'\t'gene$'\t' > temp;
      alt_search_terms=`head -n $mito_line mito_genes.txt | tail -n 1 | cut -f 5`;
      no_alts=`echo $alt_search_terms | grep -o ";" | wc -l`;
                for alt_gene in `seq 1 1 $no_alts`;
                        do alt_gene_name=`echo $alt_search_terms | cut -d ";" -f $alt_gene | sed 's/ //g'`;
                        grep "ID=gene-"$alt_gene_name";" /nesi/nobackup/uoo03398/amichael/possum_genome_master/GCF_011100635.1_mTriVul1.pri_genomic.gff | grep $'\t'gene$'\t' >> temp;
                done;
        cat temp | sort | uniq > temp_all_searches;
        if [ `wc -l temp_all_searches | awk '{ print $1 }'` -lt 1 ]
                then product_name=`grep "gene="$gene_search_term";" GCF_000001405.39_GRCh38.p13_genomic.gff | grep $'\t'CDS$'\t' | head -n 1 | sed 's/.*;product=/;product=/g' | sed 's/;protein_id=$
                grep "$product_name" /nesi/nobackup/uoo03398/michael/possum_genome_master/GCF_011100635.1_mTriVul1.pri_genomic.gff | grep ";Parent=gene-" | sed  's/.*;Parent=gene-//g' | sed 's/;.*//g' > temp_p$
                no_product_matches=`wc -l temp_product_search | awk '{ print $1 }'`;
                for product_match_line in `seq 1 1 $no_product_matches `;
                        do product_match=`head -n $product_match_line temp_product_search | tail -n 1`;
                        grep "ID=gene-"$product_match";" /nesi/nobackup/uoo03398/michael/possum_genome_master/GCF_011100635.1_mTriVul1.pri_genomic.gff | grep $'\t'gene$'\t' >> temp;
                done
                cat temp | sort | uniq > temp_all_searches;
        fi
	no_matches=`wc -l temp_all_searches | awk '{ print $1 }'`;
        echo $gene_search_term $no_matches >> results_matches.txt;
        sed "s/^/$gene_search_term\t/g" temp_all_searches >> gene_location_possum.txt;
        rm temp*;
done

## after the bash script we ran this following code to create some tables for the blast search code

module load R/4.1.0-gimkl-2020a


Type R to load R:

 

library(tidyverse)

 

results_table_tidyverse <- read_delim("results_matches.txt",delim=" ",col_names=FALSE)

 

multiple_gene_matches <- results_table_tidyverse %>% filter(X2>1)

 

zero_gene_matches <- results_table_tidyverse %>% filter(X2==0)

 

write_delim(multiple_gene_matches, "multiple_gene_matches.txt", delim=" ",col_names=FALSE)

 

write_delim(zero_gene_matches, "zero_gene_matches.txt", delim=" ",col_names=FALSE)

 

q()

# making blast databases - downloaded fna **1st april**
 # Building towards pulling out the coordinates of this gene, so we can extract the DNA and then search for it in the possum genome

grep "ID=gene-"$gene_search_term";" GCF_000001405.39_GRCh38.p13_genomic.gff | grep $'\t'gene$'\t > temp.gff

## We will grab most recent version of bedtools

module load BEDTools/2.29.2-GCC-9.2.0

## We can extract that bit of the human genome using bedtools

bedtools getfasta -fi GCF_000001405.39_GRCh38.p13_genomic.fna -bed temp.gff > temp_human.fna

module load BLAST/2.12.0-GCC-9.2.0

makeblastdb -in GCF_011100635.1_mTriVul1.pri_genomic.fna -dbtype nucl

## do for humans too

makeblastdb -in GCF_000001405.39_GRCh38.p13_genomic.fna -dbtype nucl



## next I ran this code following which incorporates two scripts (zero_match_gene_finds.sh and blast_match.R)
## when sending the job away through slurm I sent away the zero_match_gene_finds.sh slurm which within the script runs the blast_match.R code
## zero_match_gene_finds.sh script
zero_matches_total_lines=`wc -l zero_gene_matches.txt | awk '{ print $1 }'`

for zero_match_line in `seq 1 1 $zero_matches_total_lines`; 

  do gene_search_term=`head -n $zero_match_line zero_gene_matches.txt | tail -n 1 | awk '{ print $1 }'`
  
  echo $gene_search_term

  grep "gene="$gene_search_term";" GCF_000001405.39_GRCh38.p13_genomic.gff  | grep $'\t'CDS$'\t' > temp_human.gff
  
  if [ `wc -l temp_human.gff | awk '{ print $1 }'` -eq 0 ]; then
    
    echo $gene_search_term "gene_not_in_human_genome" >> blast_gene_results_matches.txt;
    
  else
  
    bedtools getfasta -fi GCF_000001405.39_GRCh38.p13_genomic.fna -bed temp_human.gff > temp_human.fna

    blastn -task blastn -db GCF_011100635.1_mTriVul1.pri_genomic.fna -query temp_human.fna -evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -outfmt 6 > temp.blast

    if [ `wc -l temp.blast | awk '{ print $1 }'` -eq 0 ]; then

      echo $gene_search_term "no_matches" >> blast_gene_results_matches.txt;

    else
    
      Rscript blast_match.R
  
      if [[ -f multiple_scaffold_matches.txt ]]; then
          mv multiple_scaffold_matches.txt $gene_search_term.multiple_scaffold_matches.txt;
          echo $gene_search_term "multiple_scaffold_matches" >> blast_gene_results_matches.txt;
      else
          if [[ -f all_blast_matches.txt ]]; then
              mv all_blast_matches.txt $gene_search_term.all_blast_matches.txt;
              echo $gene_search_term "matches_not_refined" >> blast_gene_results_matches.txt;
          else    
              found_gene_name=`bedtools intersect -a GCF_011100635.1_mTriVul1.pri_genomic.gff -b match.bed | grep $'\t'gene$'\t' | head -n 1 | sed 's/.*ID=gene-//g' | sed 's/;.*//g'`
                if [ `echo $found_gene_name | wc -w` -eq 0 ]; then
                      mv match.bed $gene_search_term.match.bed
                      echo $gene_search_term "match_not_in_gene_region" >> blast_gene_results_matches.txt;
                else
                      grep "ID=gene-"$found_gene_name";" GCF_011100635.1_mTriVul1.pri_genomic.gff | grep $'\t'gene$'\t' > temp
                      sed "s/^/$gene_search_term\t/g" temp >>  blast_match_gene_location_possum.txt
                      echo $gene_search_term $found_gene_name >> blast_gene_results_matches.txt;
                      rm match.bed
                fi      
          fi    
       fi
    fi
  fi  

  rm temp*

done

## blast_match.R script below
# Loading the library we need
library(tidyverse)

# Reading in the blast file
temp_blast <- read_delim("temp.blast",delim="\t",col_names=FALSE)

# Initialising the "while" variable
not_found_match <- TRUE

while (not_found_match) {

      # Summarising the number of matches of each possum chromosome to each human exon
      sum_matches <- temp_blast %>% group_by(X1,X2) %>% summarise(sum_match=sum(X8)) %>% pivot_wider(names_from=X2,values_from=sum_match)
      
      # Initializing the no_exons variable
      no_exons <- NULL

      # For the number of columns (i.e. possum chromosomes in sm_matches)
      for (x in 2:dim(sum_matches)[2]) {
            # Counting the number of exons that the chromosome has a match to
            no_exons <- c(no_exons,sum(!is.na(sum_matches[,x])))
      }

      # If there is a "clear winner" on number of exons matched to 
      if(sum(no_exons==max(no_exons))==1) {
            chrom_name <- names(sum_matches)[which(no_exons==max(no_exons))+1]
      } else {
      # If there isn't a clear match, look at the total sum base pair of matches      
      temp_sum <- sum_matches %>% ungroup() %>% select(which(no_exons==max(no_exons))+1) %>% colSums(na.rm=TRUE)
            # If there is a clear match based on the sum of matchs
            if(sum(temp_sum==max(temp_sum))==1) {
                  # Grab the "winner's" name 
                  chrom_name <- names(which(temp_sum==max(temp_sum)))
            } else {
                  output <- temp_blast %>% filter(X2 %in% names(temp_sum))
                  write_delim(output,"multiple_scaffold_matches.txt",delim=" ",col_names=FALSE)
                  not_found_match <- FALSE
                  break
            }
      }            

      # Whichever pathway we've taken to this point, we need to check that all the matches
      # actually belong to the same (similar) location on the scaffold
      # Match closeness is grabbing  just the matches for our winning chrom_name
      match_closeness <- temp_blast %>% filter(X2==chrom_name) %>% arrange(X9)
      # Initializing match_proximity with 1
      match_proximity <- 1
      # For each row in match_closeness starting from row 2
      for (x in 2:dim(match_closeness)[1]) {
            # If the match in row x starts no further than 75,000 bp away from the match in row x-1
            if((match_closeness$X9[x]-75000)<match_closeness$X9[x-1]) {
                  # Then give it the same value as the previous number in match_proximity
                  temp_match_proximity <- match_proximity[length(match_proximity)] 
            } else {
                  # Otherwise make it equal to the previous number in match_proximity + 1
                  temp_match_proximity <- match_proximity[length(match_proximity)]+1
            }
            # Record our temp_match_proxmity number in the match_proximity arrag
            match_proximity <- c(match_proximity,temp_match_proximity)
      }

      # If all the matches in one place
      if(length(unique(match_proximity))==1) {
         output <- match_closeness %>% 
            mutate(start=ifelse(X9>X10,X10,X9),end=ifelse(X9>X10,X9,X10)) %>% 
            select(X2,start,end) %>% unique()
         write.table(output,"match.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
         not_found_match <- FALSE
         break         
         # Otherwise, if the matches are in more than one place
      } else {
         # Recording the number of matches in each "match group"
         match_table <- as.data.frame(table(match_proximity))
         # If there is a clear winner (i.e. there isn't an even number of matches split
         # between two different locations on the chromosome
         if(length(which(match_table$Freq==max(match_table$Freq)))==1) {
            temp_blast <- temp_blast[!(do.call(paste0,temp_blast) %in% do.call(paste0,match_closeness[(match_proximity!=as.numeric(match_table$match_proximity[which(match_table$Freq==max(match_table$Freq))])),])),]
            not_found_match <- TRUE
         } else {
            write_delim(temp_blast,"all_blast_matches.txt",delim=" ",col_names=FALSE)
            not_found_match <- FALSE
            break
         }
      }
}

## the slurm script used to log the task contains the three modules that need to be loaded in order to run the code 

## created reciprocal blast script next to take the human gene name search for it in possum and then search for that back into humans to see if we get good matches completed multiple of these each using different output files from the blast script above. 
# multiple of these were created in order to use all the output files i had generated before 

# Before you run all of this, need to make human genome blast-indexed
# makeblastdb -in human_genome_name -dbtype nucl

# Name of input files the name of the input file changed we had the original one here and then blast_gene_location_possum.txt, new_gene_name_location_possum.txt, new_blast_match_gene_location_possum.txt, manual_search_list_combined.txt - i ended up making a different script separately for each of these just copying and pasting the input name differently each time
input_gene_location=gene_location_possum.txt
genome_to_pull_fasta_from=GCF_011100635.1_mTriVul1.pri_genomic
genome_to_search=GCF_000001405.39_GRCh38.p13_genomic

location_total_lines=`wc -l $input_gene_location | awk '{ print $1 }'`

for location_line in `seq 1 1 $location_total_lines`; 

  do reciprocal_blast_line=`head -n $location_line $input_gene_location | tail -n 1`
  
  human_gene_name=`echo $reciprocal_blast_line | awk '{ print $1 }'` 
  
  echo $reciprocal_blast_line | cut -f 1 -d ' ' --complement | sed 's/ /'$'\t/g' > temp_possum.gff
  
  echo $human_gene_name
  
  bedtools intersect -a $genome_to_pull_fasta_from.gff -b temp_possum.gff | grep $'\t'CDS$'\t' > temp_possum_CDS.gff
  
      if [ `wc -l temp_possum_CDS.gff | awk '{ print $1 }'` -eq 0 ]; then
      
         echo "no_CDS" $reciprocal_blast_line >> reciprocal_blast_matches.txt;
      
      else
  
        bedtools getfasta -fi $genome_to_pull_fasta_from.fna -bed temp_possum_CDS.gff > temp_possum.fna
  
        blastn -task blastn -db $genome_to_search.fna -query temp_possum.fna -evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -outfmt 6 > temp.blast

        if [ `wc -l temp.blast | awk '{ print $1 }'` -eq 0 ]; then

            echo "no_matches" $reciprocal_blast_line >> reciprocal_blast_matches.txt;

        else
    
        Rscript blast_match.R
  
          if [[ -f multiple_scaffold_matches.txt ]]; then
            mv multiple_scaffold_matches.txt $gene_search_term.rb.multiple_scaffold_matches.txt;
            echo "multiple_scaffold_matches" $reciprocal_blast_line >> reciprocal_blast_matches.txt;
         else
            if [[ -f all_blast_matches.txt ]]; then
                mv all_blast_matches.txt $gene_search_term.rb.all_blast_matches.txt;
                echo "matches_not_refined" $reciprocal_blast_line >> reciprocal_blast_matches.txt;
            else    
                found_gene_name=`bedtools intersect -a $genome_to_search.gff -b match.bed | grep $'\t'gene$'\t' | head -n 1 | sed 's/.*ID=gene-//g' | sed 's/;.*//g'`
                  if [ `echo $found_gene_name | wc -w` -eq 0 ]; then
                      mv match.bed $gene_search_term.rb.match.bed
                      echo "match_not_in_gene_region" $reciprocal_blast_line >> reciprocal_blast_matches.txt;
                  else
                      echo $found_gene_name $reciprocal_blast_line >> reciprocal_blast_matches.txt;
                      rm match.bed
                  fi      
            fi    
          fi
        fi
       
  fi  

  rm temp*

done
