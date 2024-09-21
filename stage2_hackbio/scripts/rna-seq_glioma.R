# install gplots & RColorBrewer
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

gl<-read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv", row.names = 1)
gl<-as.matrix(gl)

# Sequential Color Palette Heatmap
jpeg("glioma_heatmap_sequential.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          col = brewer.pal(9, "Blues"),  
          scale = "row",    
          key = TRUE,         
          density.info = "none", 
          trace = "none",     
          margins = c(10, 10),
	  main = "Heatmap for glioma with sequential palette")

dev.off()  

# Diverging Color Palette Heatmap
jpeg("glioma_heatmap_diverging.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          col = brewer.pal(9, "RdYlBu"),  
          scale = "row",      
          key = TRUE,         
          density.info = "none",  
          trace = "none",     
          margins = c(10, 10),
	  main = "Heatmap for glioma with diverging palette")

dev.off()  

# Clustering with sequentila & diverging colors

# sequential colors
# rows
jpeg("seq_row_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = TRUE,         
          Colv = FALSE,        
          col = brewer.pal(9, "Blues"),
          scale = "row",       
          key = TRUE,          
          density.info = "none",
	  dendrogram = "row",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by rows with seqential palette")

dev.off()

# cols
jpeg("seq_col_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = FALSE,         
          Colv = TRUE,        
          col = brewer.pal(9, "Blues"),   
          scale = "row",      
          key = TRUE,          
          density.info = "none",
	  dendrogram = "column",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by columns with sequential palette")

dev.off()

# rows & cols
jpeg("seq_both_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = TRUE,         
          Colv = TRUE,        
          col = brewer.pal(9, "Blues"),
          scale = "row",       
          key = TRUE,          
          density.info = "none",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by both rows & cols with sequential palette")

dev.off()

# diverging colors 
# rows
jpeg("div_row_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = TRUE,         
          Colv = FALSE,        
          col = brewer.pal(9, "RdYlBu"),
          scale = "row",       
          key = TRUE,          
          density.info = "none",
	  dendrogram = "row",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by rows with diverging palette")

dev.off()

# col
jpeg("div_col_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = FALSE,         
          Colv = TRUE,        
          col = brewer.pal(9, "RdYlBu"),   
          scale = "row",      
          key = TRUE,          
          density.info = "none",
	  dendrogram = "column",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by columns with diverging palette")

dev.off()

# rows & cols
jpeg("div_both_cluster.jpeg", width=8, height=8, units="in", res=300, quality=100)

heatmap.2(gl,
          Rowv = TRUE,         
          Colv = TRUE,        
          col = brewer.pal(9, "RdYlBu"),
          scale = "row",       
          key = TRUE,          
          density.info = "none",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering by both rows & cols with diverging palette")

dev.off()

# Subsetting the data: random subsetting since no metadata available
colnames(gl)

gl1<-gl[,1:5]
gl2<-gl[,6:10]

# log2FC calc
gl_fc<-rowMeans(gl1) / rowMeans(gl2)
gl_fc<-log2(gl_fc)

# p-val calc
p_val<-apply(gl, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
})

# adjusted p-val
adjusted_p_values<-p.adjust(p_val, method = "BH")

# determine cut-off : p-val = 0.5 ; fold change = 1(+), -1(-) --> results
results<-data.frame(gene=rownames(gl), LogFC=gl_fc, P_Value=p_val, Adjusted_P_Value=adjusted_p_values)

results$NegPVal<--log10(results$P_Value)

# volcano plot
plot(results$LogFC, results$NegPVal,
     pch = 20,                     
     col = ifelse(results$P_Value < 0.05 & abs(results$LogFC) > 1, "red", "black"), 
     xlab = "Log2 Fold Change", 
     ylab = "-Log10 P-Value",
     main = "Volcano Plot",
     xlim = c(-3, 3),            
     ylim = c(0, max(-log10(results$NegPVal)) + 1)) 

# draw significant threshold lines
abline(h = -log10(0.5), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

# Upregulated & downregulated genes --> DEGs
up.genes<-subset(results, Adjusted_P_Value < 0.5 & LogFC > 1)
write.csv(up.genes, "upreg_genes.csv")
down.genes<-subset(results, Adjusted_P_Value < 0.5 & LogFC < 1)
write.csv(down.genes, "downreg_genes.csv")

