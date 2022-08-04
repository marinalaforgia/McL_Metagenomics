library(microtrait)
library(piggyback)
library(tictoc)
library(parallel)
library(testthat)
detectCores()

#downloads hmms
prep.hmmmodels()

#if getting a github error need a new token
#library(gitcreds)
#gitcreds_set()

# metagenome_file <- "GL_coassembly_fixed.fa"
# microtrait_result = extract.traits(metagenome_file, save_tempfiles = T)
# saveRDS(microtrait_result, "microtrait_result.RDS")
# 
# #the result is too big
# #awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000000==0){file=sprintf("GL_coassembly_fixed_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < GL_coassembly_fixed.fa


genomes_dir="assembly_fastas"
genomes_files = list.files(genomes_dir, full.names = T, recursive = T, pattern = ".fa$")

# tictoc::tic.clearlog()
# tictoc::tic(paste0("Running microtrait for ", length(genomes_files)))
# result <- mclapply(1:length(genomes_files), function(i) {
#   r = extracttraits(genomes_files[i])
#   saveRDS(r, file = file.path(paste0(fs::path_file(genomes_files[i]), ".microtrait.rds")))
#   r}, mc.cores = 4)
# tictoc::toc(log = "TRUE")

for (i in 1:length(genomes_files)) { 
  r = extract.traits(genomes_files[i], save_tempfiles = T)
  saveRDS(r, file = file.path(paste0(fs::path_file(genomes_files[i]), ".microtrait.rds")))
  r
}


r = extract.traits(genomes_files[1], save_tempfiles = T)
