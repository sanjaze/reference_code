library(tidyverse)

get_promoters <- function(upstream_length, include_gene, cut_prevgene){
  # define promoter regions based on parameters
  # upstream length: int, bases upstream of ATG
  # include_gene: str, on of c("intron", whole", "none") considered as promoter
  # cut_prevgene: bool, TRUE: cuts off at gene in upstream length
  if (class(upstream_length) != "numeric") stop("No valid argument type for upstream_length. Must be numeric.")
  if (class(cut_prevgene) != "logical") stop("No valid argument type for cut_prevgene. Must be logical (TRUE/FALSE).")
  
  gff <- read_tsv("data_AT/TAIR10_GFF3_genes.gff", col_names = F)
  colnames(gff)<-c("chromosome","source","feature","f_start","f_end","score","strand","frame","attribute")
  
  gff <- gff %>% 
    mutate(target_id = if_else(str_detect(attribute, "(AT[:digit:]G[:digit:]+.[:digit:])"), # extract AT IDs to new column
                               str_extract(attribute, "(AT[:digit:]G[:digit:]+.[:digit:])"), 
                               str_extract(attribute, "(AT[:digit:]G[:digit:]+)"))) %>%
    mutate(note = if_else(str_detect(attribute, "Note"), str_extract(attribute, "Note=([^;]+)"), "NA")) %>%
    mutate(note = str_remove(note, "Note="))
  
  chromosomes <-pull(distinct(gff,chromosome),chromosome)
  
  new_gff <- tibble()
  
  for (chromo in chromosomes){
    genes <- gff %>%
      filter(feature == "gene" & chromosome == chromo & note == "protein_coding_gene") %>%
      arrange(f_start) # sort genes by their start on chromosome
    
    for (k in 1:nrow(genes)){ # iterate over genes on chromosome
      row <- genes[k,]
      chr1 <- row$chromosome
      
      if (include_gene == "intron"){
        exons <- gff %>% # get exons of gene
          filter(chromosome == chr1 & feature == "exon" & (f_start >= row$f_start & f_end <= row$f_end)) %>%
          arrange(chromosome, f_start)
        
        if (nrow(exons)>= 2){
          exon2_start <- exons[2,]$f_start # end of first intron on + strand
          exon2ndlast_end <- exons[nrow(exons)-1,]$f_end # end of first intron on - strand
        } else {
          exon2_start <- row$f_start
          exon2ndlast_end <- row$f_end
        }
      } else if (include_gene == "whole"){
        exon2_start <- row$f_end
        exon2ndlast_end <- row$f_start
      } else if (include_gene == "none"){
        exon2_start <- row$f_start
        exon2ndlast_end <- row$f_end
      } else {
        stop("No valid argument to include_gene. Must be either intron, whole or none.")
      }

      if (k == 1){ # first gene of chromosome
        geneprev_end <- 1
        genesubs_start <- genes[k+1,]$f_start
      } else if (k > 1 & k < nrow(genes)){
        geneprev_end <- genes[k-1,]$f_end
        genesubs_start <- genes[k+1,]$f_start
      } else {
        geneprev_end <- genes[k-1,]$f_end
        genesubs_start <- row$f_end
      }
      row_new <- cbind(row, exon2_start, exon2ndlast_end, geneprev_end, genesubs_start)
      new_gff <- bind_rows(new_gff, row_new)
    }
  }
  # define promoter regions with regards to strandedness
  if (!cut_prevgene) {
    final_gff <- new_gff %>%
      mutate(atg = if_else(strand == "+", f_start, f_end)) %>%
      mutate(prom_start = if_else(strand == "+", geneprev_end, exon2ndlast_end)) %>%
      mutate(prom_end = if_else(strand == "+", exon2_start, genesubs_start)) %>%
      mutate(upstream = if_else(strand == "+", f_start - geneprev_end, genesubs_start - f_end)) %>% # distance to upstream gene
      mutate(prom_start = if_else(strand == "+", f_start - upstream_length, prom_start)) %>% # promoter start
      mutate(prom_end = if_else(strand == "-", f_end + upstream_length, prom_end)) %>% # promoter end
      mutate(length = prom_end - prom_start)
  } else {
    final_gff <- new_gff %>%
      mutate(atg = if_else(strand == "+", f_start, f_end)) %>%
      mutate(prom_start = if_else(strand == "+", geneprev_end, exon2ndlast_end)) %>%
      mutate(prom_end = if_else(strand == "+", exon2_start, genesubs_start)) %>%
      mutate(upstream = if_else(strand == "+", f_start - geneprev_end, genesubs_start - f_end)) %>%
      # check for overlap of upstream_length with distance to upstream gene, set promoter regions acordingly:
      mutate(prom_start = if_else(upstream > upstream_length & strand == "+", f_start - upstream_length, prom_start)) %>%
      mutate(prom_end = if_else(upstream > upstream_length & strand == "-", f_end + upstream_length, prom_end)) %>%
      mutate(length = prom_end - prom_start)
  }
  
  final_gff <- final_gff %>% 
    mutate(genomic_coord = str_c(chromosome, prom_start, sep = ":")) %>%
    mutate(genomic_coord = str_c(genomic_coord, prom_end, sep = "-"))
  return(final_gff)
  
  if (FALSE){ # creates separate file for promoter regions of subset
    dir_create("out")
    write_tsv(final_gff ,"out/At_TAIR10_promoter_3kbup_intron.tsv")
    photosyn_genes <- read_tsv("data/At_PhotosynthesisGenes.tsv", col_names = T)
    photosyn_IDs <- photosyn_genes %>% select(locus)
    photosyn_prom <- merge(final_gff, photosyn_IDs, by.x = "target_id", by.y = "locus", all.x = F, all.y = T)
  
    write_tsv(photosyn_prom, "out/At_TAIR10_promoter_photosynthesis_v1.1.tsv")
  } 
}
