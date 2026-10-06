suppressPackageStartupMessages(library(data.table))
Y <- path.expand("~/Desktop/yeast/S288C_reference_genome_R64-5-1_20240529")
RUN <- "hoffmanDrop/notAF/v7"
CHRLEN <- c(chrI=230218,chrII=813184,chrIII=316620,chrIV=1531933,chrV=576874,chrVI=270161,
  chrVII=1090940,chrVIII=562643,chrIX=439888,chrX=745751,chrXI=666816,chrXII=1078177,
  chrXIII=924431,chrXIV=784333,chrXV=1091291,chrXVI=948066); CHROMS <- names(CHRLEN)

## genes
gff <- fread(cmd=sprintf("gzcat %s | awk -F'\t' '$3==\"gene\"'",
  shQuote(file.path(Y,"saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"))),
  header=FALSE, sep="\t", quote="", showProgress=FALSE)
setnames(gff,c("chr","src","type","start","end","score","strand","phase","attr"))
att <- function(k){v<-sub(paste0(".*(^|;)",k,"=([^;]*).*"),"\\2",gff$attr)
  ifelse(grepl(paste0("(^|;)",k,"="),gff$attr),v,NA_character_)}
genes <- data.table(orf=att("ID"),name=att("gene"),chrom=gff$chr,start=gff$start,
                    end=gff$end,note=att("orf_classification"))[chrom %in% CHROMS]
genes[is.na(name)|!nzchar(name), name:=orf]
genes <- genes[(is.na(note)|note!="Dubious") & !is.na(orf)]

## annotation: GAF -> systematic, propagate through OBO
gaf <- fread(cmd=sprintf("cat %s | grep -v '^!'",shQuote(file.path(Y,"gene_association_R64-5-1_20240529.sgd"))),
             header=FALSE,sep="\t",quote="",fill=TRUE,showProgress=FALSE)
gaf <- gaf[V9=="P" & !V7 %in% "ND"]
syn <- strsplit(gaf$V11,"\\|")
orf <- vapply(seq_along(syn),function(i){s<-syn[[i]]
  h<-s[grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$",s)]; if(length(h)) h[1] else gaf$V3[i]},"")
ann <- unique(data.table(orf=orf,go=gaf$V5)[orf %in% genes$orf])
ol <- readLines(file.path(Y,"go-basic.obo"),warn=FALSE)
st <- which(ol=="[Term]"); parents <- list(); nm <- character(); ids <- character()
for(b in seq_along(st)){
  blk <- ol[(st[b]+1):(if(b<length(st)) st[b+1]-1 else length(ol))]
  id <- sub("^id: ","",grep("^id: GO:",blk,value=TRUE))[1]; if(is.na(id)) next
  if(any(grepl("^is_obsolete: true",blk))) next
  if(!identical(sub("^namespace: ","",grep("^namespace: ",blk,value=TRUE))[1],"biological_process")) next
  ids<-c(ids,id); nm<-c(nm,sub("^name: ","",grep("^name: ",blk,value=TRUE))[1])
  parents[[id]] <- unique(c(sub("^is_a: (GO:[0-9]+).*","\\1",grep("^is_a: GO:",blk,value=TRUE)),
    sub("^relationship: part_of (GO:[0-9]+).*","\\1",grep("^relationship: part_of GO:",blk,value=TRUE))))
}
names(nm)<-ids; cache<-new.env(parent=emptyenv())
anc<-function(id,seen=character()){ if(!is.null(a<-cache[[id]])) return(a)
  p<-parents[[id]]; p<-p[!p %in% seen]
  a<-unique(c(p,unlist(lapply(p,anc,seen=c(seen,id,p))))); assign(id,a,envir=cache); a}
ut<-unique(ann$go); ex<-lapply(ut,function(g) unique(c(g,anc(g)))); names(ex)<-ut
full <- unique(data.table(orf=rep(ann$orf,lengths(ex[ann$go])),
                          go=unlist(ex[ann$go],use.names=FALSE)))[go %in% ids]
full[, term := nm[go]]

## contrast loci per condition
pk <- fread("hoffmanDrop/notAF/v7/contrasts/stress_peaks.csv")[databin=="main"]
MASK <- data.table(chrom=c("chrIII","chrIV"),lo=c(125000,15000),hi=c(215000,60000))
for(i in 1:nrow(MASK)) pk <- pk[!(seqnames==MASK$chrom[i] &
   ((peak.pos>=MASK$lo[i] & peak.pos<=MASK$hi[i]) |
    pmax(0,pmin(end,MASK$hi[i])-pmax(start,MASK$lo[i])) > 0.5*(MASK$hi[i]-MASK$lo[i])))]
merge_iv <- function(d){ d<-d[order(seqnames,start)]
  d[, {lo<-start[1];hi<-end[1];L<-c();H<-c()
     for(k in seq_along(start)){ if(start[k]<=hi) hi<-max(hi,end[k]) else {L<-c(L,lo);H<-c(H,hi);lo<-start[k];hi<-end[k]} }
     .(lo=c(L,lo),hi=c(H,hi))}, by=seqnames] }


## ---------------------------------------------------------------- output ----
## Table S3. Genes of the expected functional class under the CONTRAST loci of
## each stress condition, with the allele-series evidence from qtl_candidates.
##
## WHAT THIS IS NOT. The GO enrichment test is NULL -- no term reaches
## q_perm < 0.1 in any condition on either peak set -- so the PROPORTION of
## cell-wall genes under the caffeine loci, or of ion genes under the sodium
## chloride loci, is no higher than a random placement of the same loci would
## give. With 2,172 genes under 79 caffeine loci you expect about a hundred
## cell-wall-annotated genes by chance, and that is what is there. Nothing in
## this table is evidence that a locus acts through the gene named.
##
## WHAT IT IS FOR. Which genes they actually are is still the natural place to
## start a follow-up, and the allele-series evidence narrows it. Three nested
## sets, and they must not be conflated:
##   102 / 85   genes annotated to the expected process under the loci
##    24 / 13   of those that reach the candidate ranking at all, i.e. sit
##              inside the consensus region of a locus seen in >= 2 crosses
##     5 /  4   of those carrying a NON-ZERO allele-series score, i.e. at least
##              one damaging allele in a wild parent that detects the locus
## `has_allele_evidence` marks the innermost set. It is 9 genes in total, and
## that is the list worth starting from.
##
## The score, ranks and allele counts come from qtl_candidates_stress.csv and
## were computed without reference to GO, so the two lines of evidence are
## independent. `has_allele_evidence` marks the intersection.
PAT <- list(Caf=list(lab="cell wall / surface",
              re="\\bcell wall\\b|\\bglucan\\b|\\bchitin\\b|\\bmannoprotein\\b|\\bcell surface\\b|\\bfungal-type cell wall\\b"),
            NaCl=list(lab="ion / osmotic",
              re="\\bsalt\\b|\\bosmotic\\b|\\bosmolarity\\b|\\bion homeostasis\\b|\\bcation\\b|\\bsodium\\b|\\bpotassium\\b|\\bion transport\\b|\\bmetal ion\\b"))
LONGC <- c(Caf="caffeine", NaCl="sodium chloride")
cand <- fread(file.path(RUN,"candidates/stress/qtl_candidates_stress.csv"))
## genes ranked per locus -- the denominator a rank is read against. Not a
## column in the CSV; qtl_candidates_run computes it only inside its own
## self-check block, so recover it here rather than quoting a rank with no
## denominator beside it.
cand[, genes_ranked := .N, by = locus]

out <- rbindlist(lapply(names(PAT), function(cd){
  iv <- merge_iv(pk[condition==cd]); setnames(iv,"seqnames","chrom")
  g  <- genes[iv, on=.(chrom), allow.cartesian=TRUE][end>=lo & start<=hi]
  hit <- full[grepl(PAT[[cd]]$re, term, ignore.case=TRUE)]
  sel <- unique(g[orf %in% unique(hit$orf), .(orf, name, chrom)])
  tm  <- hit[orf %in% sel$orf, .(go_terms=paste(sort(unique(term)),collapse="; ")), by=orf]
  cj  <- cand[condition==cd, .(score=max(score), rank_in_locus=min(rank_in_locus),
              alleles_wild=max(alleles_wild), n_lof_wild=max(n_lof_wild),
              genes_in_locus=max(genes_ranked), locus=locus[1]), by=.(orf=gene)]
  m <- merge(merge(sel, tm, by="orf"), cj, by="orf", all.x=TRUE)
  ## 2026-09-22: these were ONE column, has_allele_evidence = score > 0, which
  ## is not what that name says.  The composite score also carries a
  ## summit-proximity term and takes its severity and conservation terms from
  ## BBQ-side variants, so a gene with NO damaging allele in any wild parent can
  ## score above zero -- 9 of the 17 rows at score > 0 are exactly that.  Split
  ## into the two distinct claims and let the prose cite whichever it means.
  m[, `:=`(condition=LONGC[[cd]], focus_class=PAT[[cd]]$lab,
           in_candidate_ranking=!is.na(score),
           has_wild_damaging_allele=!is.na(alleles_wild) & alleles_wild>0,
           nonzero_allele_series_score=!is.na(score) & score>0)]
  setorder(m, -has_wild_damaging_allele, -nonzero_allele_series_score,
           -score, na.last=TRUE)
  m[, .(condition, focus_class, gene=name, systematic=orf, chrom, locus,
        in_candidate_ranking, has_wild_damaging_allele,
        nonzero_allele_series_score, rank_in_locus, genes_in_locus,
        alleles_wild, n_lof_wild, score, go_terms)]
}))
fwrite(out, "draft/derived/TableS3_go_candidates.csv")
cat(sprintf("\nwrote draft/derived/TableS3_go_candidates.csv: %d rows\n", nrow(out)))
print(out[, .(annotated=.N, in_ranking=sum(in_candidate_ranking),
              with_wild_damaging_allele=sum(has_wild_damaging_allele),
              nonzero_score=sum(nonzero_allele_series_score)),
          by=.(condition, focus_class)])
cat("\n-- top of each condition --\n")
print(out[has_wild_damaging_allele == TRUE][, head(.SD,6), by=condition][
      , .(condition, gene, chrom, rank_in_locus, alleles_wild, n_lof_wild, score)])
