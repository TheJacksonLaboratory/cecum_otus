library(regress)
library(gap)
library(qtl2)
library(qtl2convert)
library(parallel)
library(abind)
library(ggplot2)
library(phyloseq)

metadata <- read.csv("../data/Metadata_March18_19_full.csv", row.names = 1, header=TRUE)
metadata <- metadata[metadata$Replicate != "REPLICATES",]
row.names(metadata) <- metadata$Mouse.ID
do.samples <- metadata[,c("DOB", "Sex", "Generation", "Project", "OFAdate.of.test")]
do.samples$Origin <- ifelse(do.samples$Project == "JBCC", "stool", "cecum")
GM_snps <- read.delim("../data/marker_grid_0.02cM_plus.txt", header = TRUE)
load("../data/apr_69kchr_DO1973.RData")
names(attr(apr.69kchr, "is_x_chr"))[20] = "X"

mms <- read.csv("../data/mismatch_gigamuga_39samples.csv")
probs3 = apr.69kchr
query_variants <- create_variant_query_func("../data/cc_variants.sqlite")
query_genes <- create_gene_query_func("../data/mouse_genes_mgi.sqlite")

# Compute regression
comp_reg <- function(pheno, K.overall, addcovar){
    #Assign "y" as the phenotype data extracted from the pheno object which should be a single column of data
    y = pheno
    
    t = try(regress(y~addcovar, ~K.overall, identity = TRUE), silent = T)
    
    if(!("try-error" %in% class(t))){
      #Perform regression to obtain heritability estimates
     # r = regress(y~addcovar, ~K.overall, identity = TRUE)
      r=t
      r.h2g = as.data.frame(h2G(r$sigma, r$sigma.cov, verbose = F))
      r.h2g.vp     = r.h2g$Vp
      r.h2g.vpp    = r.h2g$VVp
      r.h2g.h2g    = r.h2g$h2G
      r.h2g.varh2g = r.h2g$Varh2G
      
      r.h2g.x      = capture.output(h2G(r$sigma, r$sigma.cov, verbose = T))
      r.h2g.se     = strsplit(r.h2g.x, "SE = ")
      r.h2g.vp.se  = gsub(" ", "", r.h2g.se[[1]][2])
      r.h2g.h2g.se = gsub(" ", "", r.h2g.se[[2]][2])
      
    }else{
      r.h2g.vp = "NA"
      r.h2g.vp.se = "NA"
      r.h2g.h2g = "NA"
      r.h2g.h2g.se = "NA"
      r.h2g.varh2g = "NA"
      r.h2g.vpp = "NA"
    }
    refdf <- data.frame(r.h2g.vp, r.h2g.vp.se, r.h2g.h2g, r.h2g.h2g.se, r.h2g.varh2g, r.h2g.vpp)
    names(refdf) <- c("VP", "VP_SE", "H2G", "H2G_SE", "H2G_VAR", "VPP")
    return(refdf)
}

plot.peaks <- function(qtl.out, binob, thr = 6, model="binary"){
  # Find peaks in each coliumn of qtl.out, if peak is above threshold plot the chromosome
  for (obcol in 1:ncol(qtl.out)){
    mp <- find_peaks(qtl.out[,obcol,drop=FALSE], map)
    if (nrow(mp) > 0 ){
      for (r in 1:nrow(mp)){
        if (!is.na(mp[r,"lod"]) & (mp[r,"lod"] >= thr)){
          plot.chr(qtl.out, binob, mp[r,"chr"], obcol, model=model)
        }
      }
    }
  }
}

plot.chr <- function(qtl.out, binob, chrname, obcol, addcovar, model="binary", original=NA){
  if (is.na(original)) original = binob
  mp <- find_peaks(qtl.out[,obcol,drop=FALSE], map, drop=1.5)
  marker <- find_marker(map, chrname, pos=mp[mp$chr==chrname, "pos"])
  pname <- colnames(binob)[obcol]
  peak_Mbp <- pmap[[chrname]][marker]
  lefti <- max(peak_Mpb - 5, mp[mp$chr==chrname, "ci_lo"]) #peak_Mbp-1
  righti <- min(peak_Mbp + 5, mp[mp$chr==chrname, "ci_hi"]) #peak_Mbp+1
 # if (mp[mp$chr==chrname, "ci_hi"]-mp[mp$chr==chrname, "ci_lo"] < 4){
 #   lefti <- mp[mp$chr==chrname, "ci_lo"]
 #   righti <- mp[mp$chr==chrname, "ci_hi"]
 # }
  if (model == "binary"){
    qtl_gentcoef <- scan1coef(genoprobs[,chrname], pheno = binob[,obcol,drop=FALSE], addcovar = addcovar, model=model)
    out_snps <- scan1snps(genoprobs, pmap, binob[,obcol,drop=FALSE], addcovar=addcovar, query_func=query_variants,
                          chr=chrname, start=lefti, end=righti, keep_all_snps=TRUE, model=model)
  }else{
    qtl_gentcoef <- scan1blup(genoprobs[,chrname], pheno = binob[,obcol,drop=FALSE], addcovar = addcovar, kinship=K[[chrname]])
    out_snps <- scan1snps(genoprobs, pmap, binob[,obcol,drop=FALSE], addcovar=addcovar, kinship=K[[chrname]], query_func=query_variants,
                          chr=chrname, start=lefti, end=righti, keep_all_snps=TRUE, model=model)
  }
  plot_coefCC(qtl_gentcoef, pmap[chrname], scan1_output = qtl.out[,obcol,drop=FALSE], legend="topright", main=pname)
  genes <- query_genes(chrname, lefti, righti)
  plot_snpasso(out_snps$lod, out_snps$snpinfo, genes=genes)
  plot(qtl.out, pmap, lodcolumn=obcol, main=pname, chr=chrname)
  g <- maxmarg(genoprobs, map, chr=chrname, minprob=0.4, pos=mp[mp$chr==chrname, "pos"], return_char = TRUE)
  g2 <- maxmarg(genoprobs, map, chr=chrname, minprob=0.8, pos=mp[mp$chr==chrname, "pos"], return_char = TRUE)
  ptx<- intersect(names(g[!is.na(g)]), rownames(binob))
  px <- data.frame(g[ptx], original[ptx, obcol, drop=F], !is.na(g2[ptx]))
  colnames(px) <- c("Strain", pname, "Homo")
  print(ggplot(px, aes_string(x="Strain", y=pname, group="Strain", fill="Strain")) + geom_violin() + geom_jitter(shape=1, width = 0.1, height = 0, data=function(x){x[x$Homo==TRUE,]}) + scale_fill_manual(values=as.vector(qtl2::CCcolors), breaks=c("A","B","C","D","E","F","G","H"))+xlab(paste0("marker:", marker, " chr", chrname, ":", peak_Mbp, "Mbp")))
  ts <- top_snps(out_snps$lod, out_snps$snpinfo)
  retlist <- list("gene_names" = genes$Name, "SNP_id" = ts$snp_id[ts$lod==max(ts$lod)][1])
  return(retlist)
}


taxo <- read.delim("../data/JBCA-B-C_otu_taxa_rdpc0.5.tsv", row.names=1, sep="\t")
counts <- read.delim("../data/JBCA-B-C_otu_table_by_otutab_study-mouseID-names.txt", row.names=1)
chrorder=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
counts <- counts[,colSums(counts)>=10000]
#colnames(counts) <- sapply(sapply(colnames(counts), strsplit, "_"), "[[", 2)
for (i in 1:20){ rownames(probs3[[i]]) = sapply(sapply(rownames(probs3[[i]]), strsplit, "_"), "[[", 6)}
for (i in 1:20){ rownames(probs3[[i]])[match(as.character(mms$X...gigamuga.sample), rownames(probs3[[i]]))] <- as.character(mms$best.gbrs.sample)}
samps3 <- intersect(sapply(sapply(colnames(counts), strsplit, "_"), "[[", 2), rownames(probs3[[1]]))
#counts <- counts[,samps3]
metadata <- metadata[samps3,]
for (i in 1:20){ probs3[[i]] = probs3[[i]][samps3,,]}
genoprobs = probs3
#GM_snps <- GM_snps[abind(sapply(genoprobs, function(x) dimnames(x)[[3]])),]

#Obtain map in cM
map = map_df_to_list(map = GM_snps, pos_column="cM")
#Obtain map in bps
pmap = map_df_to_list(map = GM_snps, pos_column="pos")
apr <- genoprobs
#Specify the additive covariate model (Project is included in Flowcell)
addcovar = model.matrix(~Sex + DOB, data = do.samples)[,-1, drop = FALSE]
K = calc_kinship(probs = probs3, type = "loco", use_allele_probs = TRUE)
K.overall = calc_kinship(probs = genoprobs, type = "overall", use_allele_probs = TRUE)
# Generates a phyloseq object and agglomerate in different taxonomic levels
counts <- otu_table(counts, taxa_are_rows = TRUE)
taxo <- tax_table(as.matrix(taxo))
phys <- phyloseq(counts, taxo, metadata)
phy_gen <- tax_glom(phys, taxrank=rank_names(phys)[6])
phy_fam <- tax_glom(phys, taxrank=rank_names(phys)[5])
phy_order <- tax_glom(phys, taxrank=rank_names(phys)[4])
phy_class <- tax_glom(phys, taxrank=rank_names(phys)[3])
phy_phylum <- tax_glom(phys, taxrank=rank_names(phys)[2])
g.ord <- ordinate(phys, "NMDS", "bray", try=20)
save.image(file = "Ready_to_roll.Rdata")
