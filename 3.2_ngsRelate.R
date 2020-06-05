
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               ngsRelate               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ get file for ngsrelate

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--exclude data/processed/x_linked_variants.txt ",
              "--extract data/out/agaz/plink/agaz_plate1_3_poly_parents.prune.in ",
              #"--geno 0.1 --maf 0.01 --hwe 0.001 ",
              "--geno 0.1 --maf 0.3 --hwe 0.001 ",
              "--out data/out/agaz/plink/ngsrelate_maf ",
              "--recode vcf-iid --allow-extra-chr --debug"))

# scp data/out/agaz/plink/ngsrelate.vcf ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/
# scp data/out/agaz/plink/ngsrelate_maf.vcf ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/

#~~ read in results file from eddie

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/vcf.res data/processed/vcf.res 
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/vcf_maf.res data/processed/vcf_maf.res 

ngsrel <- fread("data/processed/vcf.res")

ggplot(ngsrel) +
  geom_point(aes(R1, R0))

ggplot(ngsrel) +
  geom_point(aes(R1, KING))

