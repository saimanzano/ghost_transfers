
library(dplyr); library(tidyr); library(pheatmap); library(gridExtra) # Imports

## Data preparing ##

# Convert data from "shiftl" object into df we can work with
results <- vector() # To store results
results2 <- vector() # To store results
for (j in (1:length(shiftl))) { # For each table (simulation)
  res <- vector() # Refresh vector
  table <- shiftl[[j]]$table # Load table
  table$donor <- gsub(".*;", "", table$donor) # Donor as last category (Phyla)
  table$sim <- j # Simulation number
  results2 <- rbind(results2, table) # Keep this results
for (i in 1:length(unique(table$iter))) { # Obtain subset
  donor1 <- table[which(table$iter == i), "donor"][1]
  donor2 <- table[which(table$iter == i), "donor"][2]
  shift <- table[which(table$iter == i), "shift"][1]
  older <- table[which(table$iter == i), "older_th_feca"][1] 
  res <- rbind(res, c(donor1, donor2, shift, older))
}
  res <- data.frame(res)
  res$simulation <- j
  results <- rbind(results, data.frame(res %>% group_by(X1,X2) %>% count(X3))) # Summarize and save
  if (j%%20 == 0) {print(j)}
}
colnames(results) <- c("donor1", "donor2", "shift", "n") # Colnames
write.table(data.frame(results), "Group_results.tsv", sep="\t") # Write results

## Figure 1: Heatmap ##

# Pair is one column ordered alphabeticallly (to cross-map info where donor1=A/donor2=B and donor1=B/donor2=A)
results$pair <- ifelse(results$donor1 < results$donor2, paste(results$donor1, results$donor2), paste(results$donor2, results$donor1))
# Group by donor pair and count shifts and no shifts
results <- data.frame(results %>% group_by(pair) %>% summarise(n_positive= sum(shift),n_negative= sum(!shift)))
# Proportion of shifts to total
results$prop <- results$n_positive/(results$n_negative + results$n_positive)
# Separate pair again
results <- results %>% separate(pair, c('donor1', 'donor2'))

# Matrix of results
mat <- spread(results[,c(1,2,5)], donor2, prop)
rownames(mat) <- mat[,1]; mat <- mat[,-1]
diag(mat) <- NA # Remove diagonal

# Heatmap
p1 <- pheatmap(mat, na_col = "white", cluster_cols = F, cluster_rows = F, border_color = NA, color = rev(hcl.colors(50, "blues3")))

## Figure 2: Donor boxplot ## 

# Again proportions but this time by shift and simulation
by_donor<- results2 %>% group_by(donor, sim) %>% summarise(n_positive= sum(shift),n_negative= sum(!shift))
by_donor$prop <- by_donor$n_positive/(by_donor$n_negative + by_donor$n_positive)
# Order by median
by_donor$donor <- factor(by_donor$donor, levels=data.frame(by_donor %>% group_by(donor) %>% summarize(median=median(prop)) %>% arrange(desc(median)))$donor)
# Boxplot
p2 <- ggplot(by_donor, aes(x=donor, y=prop, color=donor, fill=donor, alpha=0.3)) + geom_boxplot() + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1)) + xlab(element_blank()) + ylab("Proportion of shifts") 

# Arrange figures and save
cairo_pdf("Figures.pdf", width=14, height=7)
grid.arrange(grobs = list(p1[[4]],p2), ncol=2)
dev.off()


```





