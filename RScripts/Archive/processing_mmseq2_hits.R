#taking mmseqs results in & cleaning up

library(tidyverse)
library(vroom)

#ran mmseqs2 easy-search and easy-taxonomy against NR protein database
function.res.uncl <- vroom("data/mmseqs_res/functionResult_nr_GL_uncl", col_names = c("gene", "database", 
                                                                "seq_identity", "align_len",
                                                                "mismatches", "gaps",
                                                                "gene_start", "gene_end", 
                                                                "db_start", "db_end",
                                                                "evalue","bitscore"))

function.res.euk <- vroom("data/mmseqs_res/functionResult_nr_GL_euk", col_names = c("gene", "database", 
                                                                 "seq_identity", "align_len",
                                                                 "mismatches", "gaps",
                                                                 "gene_start", "gene_end", 
                                                                 "db_start", "db_end",
                                                                 "evalue","bitscore"))

#group by gene, sort by bitscore and then slice to get only top hit in database
function.res.top1.uncl <- function.res.uncl %>% # take the dataframe
  group_by(gene) %>% # group it by the grouping variable
  arrange(desc(bitscore), .by_group=TRUE) %>%
  slice(1)

function.res.top1.euk <- function.res.euk %>% # take the dataframe
  group_by(gene) %>% # group it by the grouping variable
  arrange(desc(bitscore), .by_group=TRUE) %>%
  slice(1)


write.csv(function.res.top1.uncl, "data/mmseqs_res/functional_hits_uncl.csv", quote=FALSE)
write.csv(function.res.top1.euk, "data/mmseqs_res/functional_hits_euk.csv", quote=FALSE)

#cut -f3 -d, functional_hits_euk.csv | sort | uniq > topeukhits.txt
# cut -f3 -d, functional_hits_uncl.csv | sort | uniq > topunclhits.txt

#for gene in $(cat topeukhits.txt);  do cat NRref_h | grep -w -a -m 1 $gene >> euk_resfns.txt; done
#for gene in $(cat topunclhits.txt);  do cat NRref_h | grep -w -a -m 1 $gene >> uncl_resfns.txt; done

#clean up generes.txt - add "," using sed
#remove existing ,
#sed -i '' 's/,//g' uncl_resfns.txt 
#sed -i '' 's/,//g' euk_resfns.txt 
#sed -i '' 's/ \[/,/g' uncl_resfns.txt
#sed -i '' 's/ \[/,/g' euk_resfns.txt
#sed -i '' 's/.1 /.1,/g' euk_resfns.txt
#sed -i '' 's/.1 /.1,/g' uncl_resfns.txt
#sed -i '' 's/\]/,/g' euk_resfns.txt
#sed -i '' 's/\]/,/g' uncl_resfns.txt
#mv euk_resfns.txt euk_function_headers.txt
#mv uncl_resfns.txt uncl_function_headers.txt

function_headers_E <- vroom("data/mmseqs_res/euk_function_headers.txt", col_names = c("Acession", "Function", "Taxa"))
function_headers_U <- vroom("data/mmseqs_res/uncl_function_headers.txt", col_names = c("Acession", "Function", "Taxa"))

function.res.top1.uncl <- vroom("data/mmseqs_res/functional_hits_uncl.csv", col_names = c("id","gene", "database", 
                                                                                          "seq_identity", "align_len",
                                                                                          "mismatches", "gaps",
                                                                                          "gene_start", "gene_end", 
                                                                                          "db_start", "db_end",
                                                                                          "evalue","bitscore"), skip=1)

function.res.top1.euk <-vroom("data/mmseqs_res/functional_hits_euk.csv",col_names = c("id", "gene", "database", 
                                                                                      "seq_identity", "align_len",
                                                                                      "mismatches", "gaps",
                                                                                      "gene_start", "gene_end", 
                                                                                      "db_start", "db_end",
                                                                                      "evalue","bitscore"), skip=1)

#join to get putative top hit functions
function.res.top1.function_U <- left_join(function.res.top1.uncl, function_headers_U[-c(4)], by=c("database" = "Acession"))
function.res.top1.function_E <- left_join(function.res.top1.euk, function_headers_E[-c(4)], by=c("database" = "Acession"))



#now get taxonomy last common ancestory assignments
taxonomy.lca.euk <- vroom("data/mmseqs_res/taxonomyResult_nr_GL_euk2_lca.tsv", col_names=c("gene", "databaseID", "rank", "TaxaLCA","FullLCA", "FullTopHit"), delim = "\t")
taxonomy.lca.uncl <- vroom("data/mmseqs_res/taxonomyResult_nr_GL_uncl2_lca.tsv", col_names=c("gene", "databaseID", "rank", "TaxaLCA", "FullLCA", "FullTopHit"), delim = "\t")


taxonomy.lca.euk <- taxonomy.lca.euk %>%
  mutate(FullLCA = gsub("unknown;", "", FullLCA)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms;", "", FullTopHit)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms", "unclassified", FullTopHit)) %>%
  replace(is.na(.), "unclassified") 

taxonomy.lca.uncl <- taxonomy.lca.uncl %>%
  mutate(FullLCA = gsub("unknown;", "", FullLCA)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms;", "", FullTopHit)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms", "unclassified", FullTopHit)) %>%
  replace(is.na(.), "unclassified") 
  
#join
function.res.top1.function.tax_E <- left_join(function.res.top1.function_E, taxonomy.lca.euk)
function.res.top1.function.tax_U <- left_join(function.res.top1.function_U, taxonomy.lca.uncl)

#write and assess
write.table(function.res.top1.function.tax_E, "results/function.res.top1.function.tax.euk.tsv", quote=FALSE, sep = '\t', row.names = FALSE)
write.table(function.res.top1.function.tax_U, "results/function.res.top1.function.tax.uncl.tsv", quote=FALSE, sep = '\t', row.names = FALSE)


function.summ.E <- vroom("results/function.res.top1.function.tax.euk.with.manualcuration.tsv")
function.summ.U <- vroom("results/function.res.top1.function.tax.uncl.withmanualcuration.tsv")

summary(as.factor(function.summ.E$Kingdom))
#  Bacteria    Eukaryote        Fungi        Plant Unclassified 
#53            3            6            3           17 

summary(as.factor(function.summ.U$Kingdom))
# Archaea     Bacteria    Eukaryote        Fungi     Nematode        Plant Unclassified 
# 21          412            1            1            3           38           84 

