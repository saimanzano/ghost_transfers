# Simulation of transfers
# Moisès Bernabeu
# Barcelona, May 2023

library(ggplot2)

theme_set(theme_bw())

ptr <- function(branch_space, br, ghost_birth, ghost_transfer, ghost_death) {
  p <- ggplot(brs) +
    geom_segment(aes(x = -birth, xend = -death,
                     y = reorder(1:length(brs$node), birth),
                     yend = reorder(1:length(brs$node), birth))) +
    geom_vline(xintercept = c(-feca, -leca, -(feca - leca) / 2 - leca),
               lty = 4, colour = 'steelblue') +
    xlab('Time (Mya)') +
    ylab('Branch') +
    # theme(axis.text.y = element_blank()) +
    geom_point(aes(y = as.character(br), x = -ghost_death), col = 'red') +
    geom_point(aes(y = as.character(br), x = -ghost_birth), col = 'steelblue') +
    geom_point(aes(y = as.character(br), x = -transfer), col = 'darkolivegreen')
  
  print(p)
}


sim_ghost_transfer <- function(branch_space, feca, leca) {
  # Random selection of a branch from the branch space
  br <- sample(1:dim(branch_space)[1], size = 1)
  
  # Getting the selected branch birth and death
  birth <- branch_space[br, 'birth']
  death <- branch_space[br, 'death']

  # Defining the ghosts birth, transfer and death events
  if (birth <= feca & death >= leca) {
    # If the branch is inside FECA-LECA
    ghost_birth <- runif(1, death, birth)
    transfer <- runif(1, death, ghost_birth)
    ghost_death <- runif(1, 0, transfer)
  } else if (birth >= feca & death >= leca) {
    # If the branch was born before FECA and died before LECA
    ghost_birth <- runif(1, death, birth)
    if (ghost_birth >= feca) {
      transfer <- runif(1, death, feca)
    } else {
      transfer <- runif(1, death, ghost_birth)
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
                      'transfer' = transfer)
  
  return(ghost)
}

ghost_simulator <- function(branch_space, feca, leca, N) {
  odf <- c()
  for (i in 1:N) {
    odf <- rbind(odf, sim_ghost_transfer(branch_space, feca, leca))
  }
  
  return(odf)
}

# Loading data ----
brs <- read.csv('../outputs/node_ages_phylum.tsv', sep = '\t')

# Setting boundaries
feca <- brs[which(brs$node == 'FECA-LECA'), 'birth']
leca <- brs[which(brs$node == 'FECA-LECA'), 'death']

# Filtering donor branches
brs <- brs[which(brs$birth >= leca & brs$death <= feca & brs$clades == 'Bacteria'), ]

# Simulating ghosts ----
miau <- ghost_simulator(brs, feca, leca, 1000)

# Plotting simulations
p <- ggplot(brs) +
    geom_segment(aes(x = -birth, xend = -death,
                     y = reorder(1:length(brs$node), birth),
                     yend = reorder(1:length(brs$node), birth))) +
    geom_vline(xintercept = c(-feca, -leca, -(feca - leca) / 2 - leca),
               lty = 4, colour = 'steelblue') +
    xlab('Time (Mya)') +
    ylab('Branch') +
    geom_point(data = miau, aes(y = as.character(branch), x = -ghost_death), col = 'red') +
    geom_point(data = miau, aes(y = as.character(branch), x = -ghost_birth), col = 'steelblue') +
    geom_point(data = miau, aes(y = as.character(branch), x = -transfer), col = 'darkolivegreen')
p
