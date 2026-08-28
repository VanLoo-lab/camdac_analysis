# Load metho objects
metho = qs::qread("data/out/metho/main_pure.qs")
bulk_metho = list(
  pos=metho$pos, meth=metho$meth, cov=metho$cov
)
pure_metho = list(
  pos=metho$pos, meth=metho$meth_pure, cov=metho$cov_pure
)

# Select cdif regions
chrom='1'; start=156834317; end=160189186
bpos_vec = bulk_metho$pos$chrom==chrom & bulk_metho$pos$start >= start & bulk_metho$pos$end <= end
ppos_vec = pure_metho$pos$chrom==chrom & pure_metho$pos$start >= start & pure_metho$pos$end <= end
identical(bpos_vec, ppos_vec)

# Filter objects based on vecs
bpos = bulk_metho$pos[bpos_vec,]
bmet = bulk_metho$meth[bpos_vec,]
bcov = bulk_metho$cov[bpos_vec,]
ppos = pure_metho$pos[bpos_vec,]
pmet = pure_metho$meth[bpos_vec,]
pcov = pure_metho$cov[bpos_vec,]

# Recreate metho as it was
cmetho = list(
    pos = bpos,
    meth = bmet,
    meth_pure = pmet,
    cov = bcov,
    cov_pure = pcov
)
qs::qsave(cmetho, 'data/out/metho/cdif_metho.qs')
