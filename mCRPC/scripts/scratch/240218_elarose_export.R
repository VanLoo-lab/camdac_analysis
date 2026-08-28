# Script to create folder that exports data to dropbox for Lili
library(data.table)
cln= fread('./data/pre/f4_snvmt/combined/all_valphas_cln.csv.gz')
cln = cln[patient!="DTB-255"]
cln = fwrite(cln, './data/pre/f4_snvmt/combined/all_valphas_cln.csv.gz')
av = fread('./data/pre/f4_snvmt/combined/all_valphas.csv.gz')
av = av[patient!="DTB-255"]
fwrite(av, './data/pre/f4_snvmt/combined/all_valphas.csv.gz')
ints = qs::qread('./data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_interects.qs')
str(ints,1)
ints[['DTB-255']] = NULL
qs::qsave(ints, './data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_interects.qs')
wo = fread('./data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_overlap.csv')
wo[1,]
wo = wo[sample_id!="DTB-255"]
fwrite(wo, './data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_overlap.csv')


# In bash
#dbxcli put ./data/pre/f4_snvmt/combined/all_valphas.csv.gz /VanLooP/Nana/CAMDAC_SHARE/240218_mcrpc_f123/all_valphas.csv.gz
#dbxcli put ./data/pre/f4_snvmt/combined/all_valphas_cln.csv.gz /VanLooP/Nana/CAMDAC_SHARE/240218_mcrpc_f123/all_valphas_cln.csv.gz
#dbxcli put ./data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_interects.qs /VanLooP/Nana/CAMDAC_SHARE/240218_mcrpc_f123/mcrpc_wgbs_wgs_intersect.qs
#dbxcli put ./data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_overlap.csv /VanLooP/Nana/CAMDAC_SHARE/240218_mcrpc_f123/mcrpc_wgbs_wgs_overlap.csv
