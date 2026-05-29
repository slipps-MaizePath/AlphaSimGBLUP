#!/usr/bin/env Rscript
# 01_SimulationSetup.R
# Sarah Lipps
# This is the setup script for running my multivariate simulations

library(AlphaSimR, quietly = T, verbose = F, warn.conflicts = F)
library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

print(Sys.time())
message("Setting up directories...")

#dir.create("./SimulatedData", showWarnings = FALSE)
#dir.create("./Predictions", showWarnings = FALSE)
dir.create("./Outputs", showWarnings = FALSE)
dir.create("./Outputs/Summaries", showWarnings = FALSE)

message("Directories complete. Building and optimizing the simulation grid ... ")
grid.base <- expand.grid(
  h2T1 = c(0.3, 0.5, 0.7), #maybe add 0.5 or 0.9
  h2T2 = c(0.3, 0.5, 0.7),
  h2T3 = c(0.3, 0.5, 0.7),
  rG12 = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7),
  rG13 = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7),
  rE12 = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7),
  rE13 = c(-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7)
  )

  #let's remove symmetrical duplicates
grid.new <- grid.base %>%
    mutate(det = 1-rG12^2 - rG13^2) %>%
    filter(det > 0.1) %>%
    mutate(
        T2.bundle = paste(h2T2, sprintf("%.1f", rG12), sprintf("%.1f", rE12), sep = "_"),
        T3.bundle = paste(h2T3, sprintf("%.1f", rG13), sprintf("%.1f", rE13), sep = "_")) %>%
    filter(T2.bundle >= T3.bundle) %>%
    select(-T2.bundle, -T3.bundle) %>%
    mutate(GridID = row_number())

#add missingness scenarios 
grid.new <- expand_grid(grid.new, MissingProp = c(0,0.2,0.4,0.6))


# Add seeds to the grid
N.SEEDS = 5
grid.jobs <- expand_grid(grid.new, Seed = 1:N.SEEDS) %>%
mutate(JobID = row_number()) #this part is crticical for the job array

saveRDS(grid.jobs, "./Outputs/MasterJobGrids.rds")
message(sprintf("Grid Saved! Total independent jobs: %d", nrow(grid.jobs)))

message("Simulating Master Founder Genomes...")
for(s in 1:N.SEEDS) {
  set.seed(s)
  message(paste0("Seed: ", s))
  fGenomes <- runMacs(nInd = 500, nChr = 21, inbred = FALSE, species = 'WHEAT', segSites = 8400)
  spec.SP <- SimParam$new(fGenomes)
  spec.SP$restrSegSites(overlap = F, minQtlPerChr = 100, minSnpPerChr = 500)
  spec.SP$addSnpChip(nSnpPerChr = 500)
  
  basepop <- newPop(fGenomes, simParam = spec.SP)
  fgDH <- makeDH(pop = basepop, nDH = 1, simParam = spec.SP)
  save(fgDH, spec.SP, file = paste0("./Outputs/BaseSimulation_FounderGenomes_", s, ".RData"))
}

message("Setup Complete. Ready for array submission.")
print(Sys.time())

#founderGenomes <- runMacs(nInd = 500, nChr = 21, inbred = FALSE, species = 'WHEAT', segSites = 600) #seg site = 100qtl/chr + 500 snpchip
# founderGenomes <- runMacs(nInd = 100, nChr = 5, inbred = FALSE, species = 'WHEAT', segSites = 100) #seg site = 100qtl/chr + 500 snpchip
# SP <- SimParam$new(founderGenomes)
# basePop <- newPop(founderGenomes)
# dhG <- makeDH(pop = basePop, nDH = 1)
# 
# save(founderGenomes, SP, file = "./Outputs/BaseSimulation_FounderGenomes.RData")


