
# B differential expression for manuscript

Here I am going to repeat the voom/limma analysis on 2_Transcriptome with a different significance threshold (FDR and log2(1.5)). The idea if that we have a shared B gene between A and B, we might naively expect it to be overexpressed in B+ lines by a factor of 1.5 (rather than the factor of 2 that we originally used). This will be the analysis included in the manuscript. See 2_Transcriptome for in-depth details of the analysis.

### Preprocessing

Number of genes with expression in at least one sample: 21236. 2393 genes with TPM = 0 in all samples (10.1%). NOTE: this percentage might have been misreported before, check the previous markdown files. This is due to an error in the original script -- this is the right one.

Filtering options:

```{r}
keep.exprs.group <- filterByExpr(x, group=x$samples$group,min.count=5)
keep.exprs.group[keep.exprs.group == FALSE]
x1 <- x[keep.exprs.group, keep.lib.sizes=FALSE]
dim(x1)
```
18066 genes left.

QC and MDS plots

![](manuscript/figures/suppfigD.jpeg)

### Model design

```{r}
group1=c("FB","FB","FB","MB","MB","MB","FB","FB","FB","MB","MB","MB","MB","FnoB","FnoB","FnoB","MnoB","MnoB","MnoB","FnoB","FnoB","FnoB","MnoB","MnoB","MnoB","MnoB")

design1 <- model.matrix(~0 + group1)

```

limma fit and contrast matrix

```{r}
cont.matrix1 <- makeContrasts(MB.vs.MnoB = group1MB - group1MnoB, MB.vs.FB = group1MB - group1FB, MB.vs.FnoB = group1MB - group1FnoB, FB.vs.FnoB = group1FB - group1FnoB, levels=design1)

```
Differentially expressed genes (FDR <0.05 and log2FC > 0.58)

|Genes  |MB.vs.MnoB |MB.vs.FB |MB.vs.FnoB |FB.vs.FnoB|
|-------|-----------|---------|-----------|----------|
|Down   |       147 |    2971 |      2952 |       168|
|NotSig |     23166 |   18016 |     17848 |     23216|
|Up     |       316 |    2642 |      2829 |       245|

![](misc/B_DE_for_manuscript_vennDE_comps.jpg)
![](misc/B_DE_for_manuscript_vennDE_MBvsall.jpg)

