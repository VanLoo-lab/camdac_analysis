# Calculate median coverage of mCRPC samples
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(foreach)
library(doParallel)
source('src/qc_helpers.R')

mfiles = read_lines("data/source/qc_mcrpc_methfiles.txt")
acfiles = read_lines("data/source/qc_mcrpc_acfiles.txt")
qcout = fs::dir_create("data/out/qc")
tmpdir = tempdir()

# Main - AC file counts
registerDoParallel(cores=10)
accounts = foreach(ac_file=acfiles, .combine=rbind) %dopar% {
    print(sprintf("Processing AC file %s", basename(ac_file)))
    calc_cov_boxplot(ac_file, tmpdir)
}

acresfile = path(qcout, "mcrpc_ac_cov_boxplot.tsv")
print(sprintf("Writing results to %s", acresfile))
acres = data.table::rbindlist(lapply(accounts, fread))
fwrite(acres, "data/out/qc/mcrpc_ac_cov_boxplot.tsv", sep="\t")
for(i in accounts){ fs::file_delete(i) }

# Main - Meth file counts
registerDoParallel(cores=10)
methcounts = foreach(meth_file=mfiles, .combine=rbind) %dopar% {
    print(sprintf("Processing meth file %s", basename(meth_file)))
    calc_n_meth_pass(meth_file, tmpdir)
}

methresfile = path(qcout, "mcrpc_meth_n_meth_pass.tsv")
print(sprintf("Writing results to %s", methresfile))
methres = data.table::rbindlist(lapply(methcounts, fread))
fwrite(methres, methresfile, sep="\t")
for(i in methcounts){ fs::file_delete(i) }
