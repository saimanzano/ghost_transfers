# Simulation of transfers
# Moisès Bernabeu
# Barcelona, May 2023

library(ggplot2)
library(tidyr)
library(ggpubr)

theme_set(theme_bw())

ptr <- function(branch_space, sim_res = NULL, sim_res_2 = NULL, iter = NULL, sorting = 'birth') {
  if (!is.null(sim_res_2)) {
    if (!is.null(iter)) {
      sim_res <- rbind(sim_res[iter, ], sim_res_2[iter, ])
    } else {
      sim_res <- rbind(sim_res, sim_res_2)
    }
  }
  
  p <- ggplot(branch_space) +
    geom_segment(aes(x = -birth, xend = -death,
                     y = reorder(1:length(branch_space$node), get(sorting)),
                     yend = reorder(1:length(branch_space$node), get(sorting)))) +
    geom_vline(xintercept = c(-feca, -leca),
               lty = 4, colour = 'steelblue') +
    xlab('Time (Mya)') +
    ylab('Branch')
  
    if (!is.null(sim_res)) {
      p <- p + geom_point(data = sim_res, aes(y = as.character(branch), x = -ghost_death), col = 'red', alpha = 0.4) +
        geom_point(data = sim_res, aes(y = as.character(branch), x = -ghost_birth), col = 'steelblue', alpha = 0.4) +
        geom_point(data = sim_res, aes(y = as.character(branch), x = -transfer), col = 'darkolivegreen', alpha = 0.4)
    }
  
  print(p)
}

ghost_transfer <- function(branch_space, feca, leca, weights = NULL) {
  # Random selection of a branch from the branch space
  if (is.null(weights)) {
    br <- sample(1:dim(branch_space)[1], size = 1)
  } else {
    br <- sample(1:dim(branch_space)[1], size = 1, prob = weights)
  }
  donor <- branch_space[br, 'phylums']

  # Getting the selected branch birth and death
  birth <- branch_space[br, 'birth']
  death <- branch_space[br, 'death']

  # Defining the ghosts birth, transfer and death events
  if (birth <= feca & death >= leca) {
    # If the branch is inside FECA-LECA
    ghost_birth <- runif(1, death, birth)
    transfer <- runif(1, leca, ghost_birth)
    ghost_death <- runif(1, 0, transfer)
  } else if (birth >= feca & death >= leca) {
    # If the branch was born before FECA and died before LECA
    ghost_birth <- runif(1, death, birth)
    if (ghost_birth >= feca) {
      transfer <- runif(1, leca, feca)
    } else {
      transfer <- runif(1, leca, ghost_birth)
    }
    ghost_death <- runif(1, 0, transfer)
  } else if (birth <= feca & death <= leca) {
    # If the branch was born after FECA and died after LECA
    ghost_birth <- runif(1, leca, birth)
    transfer <- runif(1, leca, ghost_birth)
    ghost_death <- runif(1, 0, transfer)
  } else {
    # If the branch was born after FECA and died after LECA
    ghost_birth <- runif(1, leca, birth)
    transfer <- runif(1, leca, feca)
    ghost_death <- runif(1, 0, transfer)
  }

  # Returning the information
  ghost <- data.frame('branch' = br,
                      'sister_birth' = birth, 'sister_death' = death,
                      'ghost_birth' = ghost_birth, 'ghost_death' = ghost_death,
                      'transfer' = transfer, 'donor' = donor)
  
  return(ghost)
}

ghost_simulator <- function(branch_space, feca, leca, N, weights = NULL) {
  odf <- c()
  for (i in 1:N) {
    odf <- rbind(odf, ghost_transfer(branch_space, feca, leca, weights))
  }
  return(odf)
}

get_shifts <- function(sim_1, sim_2, feca) {
  if (dim(sim_1)[1] != dim(sim_2)[1]) {
    stop('Both simulations have not the same number of iteration, they are not comparable.')
  }
  
  out <- c()
  for (iter in 1:1000) {
    pair <- rbind(sim_1[iter, ], sim_2[iter, ])
    mintrans <- which.min(pair$transfer)
    minbirth <- which.min(pair$ghost_birth)
    pair[mintrans, 'transfer_comparison'] <- 'early'
    pair[-mintrans, 'transfer_comparison'] <- 'late'
    pair[minbirth, 'conclusion'] <- 'early'
    pair[minbirth, 'shift'] <- pair[minbirth, 'transfer_comparison'] != pair[minbirth, 'conclusion']
    pair[-minbirth, 'conclusion'] <- 'late'
    pair[-minbirth, 'shift'] <- pair[-minbirth, 'transfer_comparison'] != pair[-minbirth, 'conclusion']
    pair[minbirth, 'older_th_feca'] <- pair[minbirth, 'ghost_birth'] > feca
    pair[-minbirth, 'older_th_feca'] <- pair[-minbirth, 'ghost_birth'] > feca
    pair[, 'iter'] <- iter

    row.names(pair) <- NULL
    out <- rbind(out, pair)
    # print(iter)
  }
  
  summ <- table(out$shift) / dim(out)[1] * 100
  names(summ) <- c('No shift', 'Shift')
  
  feca_summ <- table(out$older_th_feca) / dim(out)[1] * 100
  names(feca_summ) <- c('More recent than FECA', 'Older than FECA')

  shifts_id <- which(out$shift == TRUE)
  shift_and_older <- table(out[shifts_id, 'older_th_feca']) / length(shifts_id) * 100
  names(shift_and_older) <- c('Shift more recent than FECA', 'Shift older than FECA')

  return(list('table' = out, 'shif_summary' = summ, 'out_of_feca_summary' = feca_summ,
              'shift_and_older' = shift_and_older))
}

# Loading data ----
brs <- read.csv('../outputs/node_ages_phylum.tsv', sep = '\t')

# Setting boundaries
feca <- brs[which(brs$node == 'FECA-LECA'), 'birth']
leca <- brs[which(brs$node == 'FECA-LECA'), 'death']

# Filtering donor branches
brs <- brs[which(brs$birth >= leca & brs$death <= feca & brs$clades == 'Bacteria'), ]
brs$node <- factor(brs$node)

# Choosing weights for the branch selection
brs[which(brs$death <= leca & brs$birth >= feca), 'feca_leca_prop'] <- 1
left_censored <- which(brs$death <= leca & brs$birth <= feca)
right_censored <- which(brs$death >= leca & brs$birth >= feca)
brs[left_censored, 'feca_leca_prop'] <- (brs[left_censored, 'birth'] - leca) / (feca - leca)
brs[right_censored, 'feca_leca_prop'] <- (feca - brs[right_censored, 'death']) / (feca - leca)
inside_feca_leca <- which(brs$death >= leca & brs$birth <= feca)
brs[inside_feca_leca, 'feca_leca_prop'] <- (brs[inside_feca_leca, 'birth'] - brs[inside_feca_leca, 'death']) / (feca - leca)

prefeca <- which(brs$birth >= feca)
brs[prefeca, 'prefeca_prop'] <-(brs[prefeca, 'birth'] - feca) / (max(brs[prefeca, 'birth']) - feca)
brs[which(brs$birth < feca), 'prefeca_prop'] <- 0
brs[, 'score'] <- (brs[, 'prefeca_prop'] + brs[, 'feca_leca_prop']) / max(brs[, 'prefeca_prop'] + brs[, 'feca_leca_prop'])
brs[, 'weight'] <- brs[, 'score'] / sum(brs[, 'score'])

a <- ptr(brs)
b <- ptr(brs, sorting = 'weight') +
  geom_point(aes(x = -death, y = reorder(1:length(brs$node), weight), colour = weight)) +
  scale_colour_gradient(high = '#132B43', low = '#56B1F7') +
  theme(legend.text = element_text(angle = 45, hjust = 1))

ggarrange(a, b, common.legend = TRUE)

run_1 <- ghost_simulator(brs, feca, leca, 1000, brs[, 'weight'])
run_2 <- ghost_simulator(brs, feca, leca, 1000, brs[, 'weight'])

ptr(brs, run_1, run_2)

# Conclusion shift test
shifts <- get_shifts(run_1, run_2, feca)
shifts

shiftl <- c()
tm0 <- Sys.time()
for (i in 1:1000) {
  shiftl[[i]] <- get_shifts(ghost_simulator(brs, feca, leca, 1000, brs[, 'weight']),
                            ghost_simulator(brs, feca, leca, 1000, brs[, 'weight']),
                            feca)
  print(i)
}
tm1 <- Sys.time()
tm1 - tm0

shift <- gather(data.frame(do.call(rbind, lapply(shiftl, function(x) {x$shif_summary})), check.names = FALSE))
shift_older <- gather(data.frame(do.call(rbind, lapply(shiftl, function(x) {x$shift_and_older})), check.names = FALSE))

shift_older$key_2 <- factor(shift_older$key, labels = c(FALSE, TRUE))

a <- ggplot(shift, aes(x = key, y = value, colour = key, fill = key)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Shift') +
  ylim(10, 90)

b <- ggplot(shift_older, aes(x = key_2, y = value, colour = key, fill = key)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Shift older than FECA') +
  ylim(10, 90)

pdf('../outputs/shift_older_proportions_weighted.pdf', width = 5, height = 4)
ggarrange(a, b, align = 'hv', legend = FALSE)
dev.off()

# Actinobacteria - alphaproteobacteria
actalphbrs <- brs[which(brs$class == 'Alphaproteobacteria' | brs$phylums == 'Actinobacteria'), ]
ptr(actalphbrs) +
  geom_point(aes(x = -death, y = reorder(1:length(actalphbrs$node), birth), colour = class))

actal_run_1 <- ghost_simulator(brs[which(brs$class == 'Alphaproteobacteria'), ], feca, leca, 1000)
actal_run_2 <- ghost_simulator(brs[which(brs$phylums == 'Actinobacteria'), ], feca, leca, 1000)

actal_shifts <- get_shifts(actal_run_1, actal_run_2, feca)

actal_shifts
