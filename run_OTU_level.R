library(regress)
library(gap)
library(qtl2)
library(qtl2convert)
library(parallel)
library(abind)
library(ggplot2)
library(phyloseq)
alldata = data.frame()
allh2g = data.frame()
# treat all numbers <10 as zeros and compute relative abundance. Use it in normal mode
if (file.exists("Ready_to_roll_OTU.Rdata") && (file.mtime("Ready_to_roll_OTU.Rdata") > file.mtime("Ready_to_roll.Rdata"))) {
    load("Ready_to_roll_OTU.Rdata")
}else{
    load("Ready_to_roll.Rdata")
    relgenus = as.data.frame(base::t(otu_table(phys)))
    relgenus[relgenus<10] <- NA
    relgenus <- relgenus[, colSums(!is.na(relgenus))>=10]
    origgenus <- as.data.frame(base::t(otu_table(phys)))
    origgenus <- origgenus[,colSums(!is.na(relgenus))>=10]
    cori <- origgenus[substr(row.names(origgenus),1,4)!="JBCC",]
    sori <- origgenus[substr(row.names(origgenus),1,4)=="JBCC",]
    row.names(cori) <- sapply(sapply(row.names(cori), strsplit, "_"), "[[", 2)
    row.names(sori) <- sapply(sapply(row.names(sori), strsplit, "_"), "[[", 2)
    rname <- row.names(relgenus)
    relgenus <- as.data.frame(apply(relgenus, 2, function(x) as.numeric(x/sum(x, na.rm=T))))
    row.names(relgenus) <- rname
    crel <- relgenus[substr(row.names(relgenus),1,4)!="JBCC",]
    srel <- relgenus[substr(row.names(relgenus),1,4)=="JBCC",]
    row.names(crel) <- sapply(sapply(row.names(crel), strsplit, "_"), "[[", 2)
    row.names(srel) <- sapply(sapply(row.names(srel), strsplit, "_"), "[[", 2)
#    qtl.genus2.c = scan1(genoprobs = genoprobs, pheno = crel, kinship = K, addcovar = addcovar, model="normal", cores=32)
#    qtl.genus2.s = scan1(genoprobs = genoprobs, pheno = srel, kinship = K, addcovar = addcovar, model="normal", cores=32)
    sqgenus = as.data.frame(base::t(sqrt(otu_table(phys))))
#    colnames(sqgenus) <- paste0('g__', gsub("[ /]", "_", tax_table(phy_gen)[,6]))
    sqgenus <- sqgenus[,colnames(relgenus)]
    csq <- sqgenus[substr(row.names(sqgenus),1,4)!="JBCC",]
    ssq <- sqgenus[substr(row.names(sqgenus),1,4)=="JBCC",]
    row.names(csq) <- sapply(sapply(row.names(csq), strsplit, "_"), "[[", 2)
    row.names(ssq) <- sapply(sapply(row.names(ssq), strsplit, "_"), "[[", 2)
    addcovar.c <- addcovar[intersect(rownames(csq), rownames(addcovar)),]
    addcovar.c <- cbind(addcovar.c[,colSums(addcovar.c>0)>0], data.frame(depth=rowSums(csq[rownames(addcovar.c),])))
    addcovar.s <- addcovar[intersect(rownames(ssq), rownames(addcovar)),]
    addcovar.s <- cbind(addcovar.s[,colSums(addcovar.s>0)>0], data.frame(depth=rowSums(ssq[rownames(addcovar.s),])))
    print(head(addcovar.c))
    print(head(addcovar.s))
    print(head(csq))
    print(head(ssq))
    print(dim(csq))
    print(dim(ssq))
    qtl.genus.sq.c = scan1(genoprobs = genoprobs, pheno = csq, addcovar = addcovar.c, model="normal", kinship=K, cores=1)
    print("computed C")
    qtl.genus.sq.s = scan1(genoprobs = genoprobs, pheno = ssq, addcovar = addcovar.s, model="normal", kinship=K, cores=1)
    save.image("Ready_to_roll_OTU.Rdata")
}
args = commandArgs(trailingOnly=TRUE)
ntest = args[1]
rownum = args[2]
print(ntest)
print(rownum)
pdf(paste0("OTU_QTLs/OTU_level_analysis_all_trans_p1000_",ntest,"_", rownum, "_withJBCC.pdf"))

if (ntest=="stool"){
 #   qtl.genus2 <- qtl.genus2.s
    qtl.genus.sq <- qtl.genus.sq.s
    sqgenus <- ssq
    relgenus <- srel
    origgenus <- sori
    addcovar <- addcovar.s
  }else if(ntest=="cecum"){
  #  qtl.genus2 <- qtl.genus2.c
    qtl.genus.sq <- qtl.genus.sq.c
    sqgenus <- csq
    relgenus <- crel
    origgenus <- cori
    addcovar <- addcovar.c
  }else{
    print ("Please specify 'stool' or 'cecum'")
    quit()
  }
print(addcovar)

for (i in c(as.numeric(rownum))){
  print(i)
  print(colnames(sqgenus)[i])
  print ("Normal model, sqrt transform")
  rn = intersect(intersect(row.names(K.overall), row.names(sqgenus)), row.names(addcovar))
#    ag <- comp_reg(sqgenus[rn,i], K.overall[rn, rn], addcovar[rn,])
  aag <- est_herit(sqgenus[rn,i, drop=FALSE], K.overall[rn, rn], addcovar[rn,])
  ag <- data.frame(heritability=aag, log10likelihood=attr(aag, "log10lik"))
  ag$method = "sqrt"
  ag$genus = colnames(sqgenus)[i]
  ag$origin <- ntest
  allh2g <- rbind(allh2g, ag)
  if (!is.na(max(qtl.genus.sq[,i])) && !is.infinite(max(qtl.genus.sq[,i]))){
#    mp <- find_peaks(qtl.genus.sq[,i,drop=FALSE], map, threshold=5)
    mpp <- find_peaks(qtl.genus.sq[,i,drop=FALSE], pmap, threshold=5, drop=1.5)
    if (nrow(mpp)>0){
      mpp$pvalue <- 1
      mpp$method <- "sqrt"
      mpp$origin <- ntest
      mpp$marker <- sapply(rownames(mpp), function(x) find_marker(pmap, mpp[x, "chr"], pos=mpp[x, "pos"]))
      mpp$nonzero <- colSums(sqgenus>0)[i]
    }
    adds = ""
    if ((max(mpp$lod) > 5) && (colSums(sqgenus>0)[i]>0)){
      ap <- scan1perm(apr, sqgenus[,i,drop=FALSE], model="normal", kinship=K, n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(map), addcovar=addcovar, cores=1)
      print(summary(ap, alpha=c(0.01, 0.05, 0.1, 0.2)))
      adds = sprintf(", p-value threshold: %g",summary(ap, 0.05)$A["0.05",])
      hist(ap$A, main=paste0(colnames(sqgenus)[i], " normal model - sqrt, permutations"), breaks=50)
      if (nrow(mpp)>0){
        mpp$pvalue = sapply(rownames(mpp), function(x) if (mpp[x, "chr"]=="X"){sum(ap$X >= mpp[x,"lod"])/length(ap$X)}else{sum(ap$A >= mpp[x, "lod"])/length(ap$A)})
      }
    }
    hist(sqgenus[,i], breaks = 50, main=paste("Sqrt of read counts-", colnames(sqgenus)[i], adds, sep=""))
    plot(qtl.genus.sq, map, lodcolumn=i, main=paste("Sqrt of read counts-", colnames(sqgenus)[i], adds, sep=""), chr=chrorder)
    for (r in 1:nrow(mpp)){
      if (mpp[r, "lod"] >= 5){
        rl <- plot.chr(qtl.genus.sq, sqgenus, mpp[r, "chr"], i, addcovar, model="normal", original=origgenus)
        mpp[r, "SNP_id"] <- rl$SNP_id
        mpp[r, "gene_names"] <- paste(rl$gene_names, collapse=";")
      }
    }
    if (nrow(mpp)>0){
      alldata <- rbind(alldata, mpp)
    }
  }
}
write.csv(alldata, paste0("OTU_QTLs/all_data_OTU_level_",ntest,"_",rownum, ".csv"))
write.csv(allh2g, paste0("OTU_QTLs/all_h2g_calc_OTU_level_",ntest,"_",rownum, ".csv"))

dev.off()
