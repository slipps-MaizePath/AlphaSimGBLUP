#!/usr/bin/env Rscript
# 02_RunPipeline.R
# Sarah Lipps
# 20260518
# This script runs the simulation functions using the trait grid and founder genomes from 01_SimulationSetup.R

library(AlphaSimR)
library(asreml)
library(tidyverse)
source("./MasterUtilityFunctions.R")


# this is designed to work with a SLURM job on the biocluster but will still function on a local machine
CHUNK_SIZE <- 1000
#CHUNK_SIZE <- 1
ARRAY_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(ARRAY_ID)) {
    warning("No SLURM_ARRAY_TASK_ID found. Defaulting to Job 1 for local testing.")
  ARRAY_ID = 1
}

# load objects generated in 01_SimulationSetup
JobGrid <- readRDS("./Outputs/MasterJobGrids.rds") #grid of parameters
#load("./Outputs/BaseSimulation_FounderGenomes.RData") # 'founderGenomes' and 'SP'

START.IDX <- (ARRAY_ID - 1) * CHUNK_SIZE + 1
END.IDX <- min(ARRAY_ID * CHUNK_SIZE, nrow(JobGrid))

if(START.IDX > nrow(JobGrid)) stop("Invalid Task ID. Halting Pipeline.")
p_chunk <- JobGrid[START.IDX:END.IDX, ]
message(paste("=== STARTING CHUNK: ", ARRAY_ID, " | PROCESSING ROWS: ", START.IDX, " TO ", END.IDX))
chunk.summaries <- list()

for (row_i in 1:nrow(p_chunk)) {
  p <- p_chunk[row_i, ]
  combo.tag <- sprintf("GRID%05d_Seed%02d", p$JobID, p$Seed)
  set.seed(20260518 + p$JobID + p$Seed)
  message(paste0(" -> Running Simulation: ", combo.tag))
  
  # Load the Foundation
  #foundation.path <- paste0("./Outputs/Foundation_Seed_", p$Seed, ".RData")
  foundation.path <- paste0("./Outputs/BaseSimulation_FounderGenomes_", p$Seed, ".RData")
  load(foundation.path)
  SP <- spec.SP
  #Simulate Traits and apply missingness
  sim.list <- SimulateTraitsMarkers(
    SP = SP, POP = fgDH, 
    rG12 = p$rG12, rG13 = p$rG13, rG23 = 0.0,
    rE12 = p$rE12, rE13 = p$rE13, rE23 = 0.1,
    h2T1 = p$h2T1, h2T2 = p$h2T2, h2T3 = p$h2T3,
    MISSING.PROP = p$MissingProp, # Ensure this column exists in your MasterJobGrid
    nreps = 3
  )
  
  my.pop = sim.list$pop
  my.pheno = sim.list$pheno.df
  
  #store the True Breeding values for accuracy calculation
  T.BV <- data.frame(
    id = as.factor(my.pop@id),
    Yield.TBV = bv(my.pop)[,1]
  )
  
  SNP.MATRIX <- pullSnpGeno(my.pop, simParam = spec.SP)
  rownames(SNP.MATRIX) <- my.pop@id
  Gobj <- MakeGRM(SNP.MATRIX = SNP.MATRIX)
  Ginv <- Gobj$G.inv
  GRM <- Gobj$GRM
  
 # RUN MODELS
  out.uniF <- m.asr(PHENO.DAT = my.pheno, TRAITS.LIST = "Yield", OPTION = "uni", USE.GRM = F)
  out.uniT <- m.asr(GRM= Ginv, PHENO.DAT = my.pheno, TRAITS.LIST = "Yield", OPTION = "uni", 
		    USE.GRM = T)
  out.bivF <- m.asr(PHENO.DAT = my.pheno, TRAITS.LIST = list(c("Yield", "Standard2o"), 
		    c("Yield", "Latent2o")), OPTION = "biv", USE.GRM = F, GRM = NULL)
  out.bivT <- m.asr(GRM = Ginv, PHENO.DAT = my.pheno, 
		    TRAITS.LIST = list(c("Yield","Standard2o"),c("Yield", "Latent2o")), 
                    OPTION = "biv", USE.GRM = T)
  out.triF <- m.asr(PHENO.DAT = my.pheno, 
		    TRAITS.LIST = list(c("Yield", "Standard2o", "Latent2o")), 
                    OPTION = "tri", USE.GRM = F)
  out.triT <- m.asr(GRM = Ginv, PHENO.DAT = my.pheno, 
		    TRAITS.LIST = list(c("Yield", "Standard2o", "Latent2o")), 
                    OPTION = "tri", USE.GRM = T)
  
  #get statistics all in one go
  uniF <- ExtractUni(out.uniF, USE.GRM=F, MODEL.NAME="Uni_Standard", SIM.NUM=p$JobID)
  uniT <- ExtractUni(out.uniT, USE.GRM=T, GRM=GRM, MODEL.NAME="Uni_GBLUP", SIM.NUM=p$JobID)
  bivT <- ExtractMulti(out.bivT, USE.GRM=T, GRM=GRM, MODEL.NAME="Biv_GBLUP", SIM.NUM=p$JobID)
  bivF <- ExtractMulti(out.bivF, USE.GRM=F, MODEL.NAME="Biv_Standard", SIM.NUM=p$JobID)
  triT <- ExtractMulti(out.triT, USE.GRM=T, GRM=GRM, MODEL.NAME="Tri_GBLUP", SIM.NUM=p$JobID)
  triF <- ExtractMulti(out.triF, USE.GRM=F, MODEL.NAME="Tri_Standard", SIM.NUM=p$JobID)

  summaries <- do.call(rbind, list(uniF$summary.df, uniT$summary.df, bivT$summary.df,
                                  bivF$summary.df, triT$summary.df, triF$summary.df))

  blups <- do.call(rbind, list(uniF$blups, uniT$blups, bivT$blups,
			       bivF$blups, triT$blups, triF$blups))  

  emp.accy <- blups %>%
	filter(trait %in% c("Yield", "trait_Yield")) %>%
	left_join(T.BV, by = "id") %>%
	group_by(ModelType, fold) %>%
	summarise(
	  Emp.Accy.Yield = cor(Yield.TBV, predicted.value, use = "complete.obs"),.groups="drop")

  final.df <- summaries %>%
	left_join(emp.accy, by = c("ModelName"="ModelType"))
  
  chunk.summaries[[row_i]] <- final.df
}

final.chunk.output <- bind_rows(chunk.summaries)
write.csv(final.chunk.output, paste0("./Outputs/Summaries/Summary_Chunk_", ARRAY_ID, ".csv"), row.names = F)

message(paste0("=== SUCCESS: COMPLETED CHUNK ", ARRAY_ID, " ==="))
# 
# # get parameters for this specific job
# p <- JobGrid[JobGrid$JobID == TASK_ID, ]
# if(nrow(p) == 0) {
#     stop("Invalid Task ID. Halting pipeline.")
# }
# 
# #i <- p$GridID
# #seed <- p$Seed
# combo.tag <- sprintf("GRID%05d_Seed%02d", p$GridID, p$Seed)
# set.seed(20260518 + p$GridID + p$Seed)
# 
# message(paste0("=== STARTING JOB: ", TASK_ID, " | TAG: ", combo.tag, " ==="))
#     
#     ###############################################
#     # SIMULATE TRAITS AND EXTRACT BREEDING VALUES #
#     ############################################### 
#     # simulate
#     sim.list <- SimulateTraitsMarkers(
#       SP = SP, GENOMES = founderGenomes, 
#       rG12 = p$rG12, rG13 = p$rG13, rG23 = 0.0,
#       rE12 = p$rE12, rE13 = p$rE13, rE23 = 0.1,
#       h2T1 = p$h2T1, h2T2 = p$h2T2, h2T3 = p$h2T3,
#       nreps = 3
#     )
#     
#     my.pop <- sim.list[[1]]
#     my.pheno <- sim.list[[2]]
#     
#     # get the breeding values
#     T.BV <- data.frame(
#       id = as.factor(my.pop@id),
#       Yield.TBV = bv(my.pop)[, 1],
#       Standard2o.TBV = bv(my.pop)[, 2],
#       Latent2o.TBV = bv(my.pop)[, 3],
#       Seed = p$Seed,
#       ComboTag = combo.tag
#     )
#     #save the simulated data
#     #saveRDS(T.BV, file = paste0("./SimulatedData/Sim_", combo.tag, "_TBV.rds"))
#     #saveRDS(my.pheno, file = paste0("./SimulatedData/Sim_", combo.tag, "_Phenotypes.rds"))
# 
#     ################
#     # B CREATE GRM #
#     ################
#     SNP.MATRIX <- pullSnpGeno(my.pop)
#     rownames(SNP.MATRIX) <- my.pop@id
#     
#     Glist <- MakeGRM(SNP.Matrix = SNP.MATRIX)
#     #G.inv <- Glist$G.inv
#     #GRM <- Glist$GRM
# 
#     ################################
#     # C RUN ASREML & EXTRACT BLUPS #
#     ################################  
#     
#     #univariate
#     uni.traits <- c("Yield")
#     uni.grmT <- m.asr(PHENO.DAT = my.pheno, OPTION = "uni", TRAITS.LIST = uni.traits,
#                       GRM = Glist$G.inv, USE.GRM = T)
#     out.uniT <- ExtractMaster(uni.grmT, USE.GRM = T, GRM = Glist$GRM, MODEL.NAME = "Uni_G",
#                               SIM.NUM = combo.tag)
#     uni.grmF <- m.asr(PHENO.DAT = my.pheno, OPTION = "uni", TRAITS.LIST = uni.traits)
#     out.uniF <- ExtractMaster(uni.grmF, USE.GRM = F, GRM = NULL, MODEL.NAME = "Uni_A",
#                               SIM.NUM = combo.tag)
#     
#     #bivariate (Yield + Standard 2o and Yield + Latent 2o)
#     bi.traits <- list(c("Yield", "Standard2o"), c("Yield", "Latent2o"))
#     biv.grmT <- m.asr(PHENO.DAT = my.pheno, OPTION = "biv", TRAITS.LIST = bi.traits,
#                       GRM = Glist$G.inv, USE.GRM = T)
#     out.bivT <- ExtractMaster(biv.grmT, USE.GRM = T, GRM = Glist$GRM, MODEL.NAME = "Biv_G",
#                               SIM.NUM = combo.tag)
#     biv.grmF <- m.asr(PHENO.DAT = my.pheno, OPTION = "biv", TRAITS.LIST = bi.traits)
#     out.bivF <- ExtractMaster(biv.grmF, USE.GRM = F, GRM = NULL, MODEL.NAME = "Biv_A",
#                               SIM.NUM = combo.tag)
#     
#     #trivariate (Yield + Standard2o + Latent 2o)
#     tri.traits <- list(c("Yield", "Standard2o", "Latent2o"))
#     tri.grmT <- m.asr(PHENO.DAT = my.pheno, OPTION = "tri", TRAITS.LIST = tri.traits,
#                       GRM = Glist$G.inv, USE.GRM = T)
#     out.triT <- ExtractMaster(tri.grmT, USE.GRM = T, GRM = Glist$GRM, MODEL.NAME = "Tri_G",
#                               SIM.NUM = combo.tag)
#     tri.grmF <- m.asr(PHENO.DAT = my.pheno, OPTION = "tri", TRAITS.LIST = tri.traits)
#     out.triF <- ExtractMaster(tri.grmF, USE.GRM = F, GRM = NULL, MODEL.NAME = "Tri_A",
#                               SIM.NUM = combo.tag)
# 
#    #combine blups from all models
#     currrent.sim.blups <- bind_rows(out.uniF$blups, out.uniT$blups, out.bivF$blups,
#                                     out.bivT$blups, out.triF$blups, out.triT$blups)
#     
#     #combine summaries
#     current.sim.summaries <- bind_rows(out.uniF$summary.df, out.uniT$summary.df, out.bivF$summary.df,
#                                        out.bivT$summary.df, out.triF$summary.df, out.triT$summary.df) %>%
#       dplyr::mutate(GridID = p$GridID, Seed = p$Seed) %>%
#       dplyr::bind_cols(p %>% select(-any_of(c("JobID", "Seed", "GridID")))) #this attaches the true grid parameters to each row
#     
#     #calculating the empirical accuracy by correlating predicted bvs with true bvs
#     emp.accy <- currrent.sim.blups %>%
#       filter(trait == "Yield" | trait == "trait_Yield") %>% #only looking at yield blups
#       left_join(T.BV, by = "id") %>%
#       group_by(ModelType, fold) %>%
#       summarise(Empirical.Accuracy.Yield = cor(Yield.TBV, predicted.value, use = "complete.obs"),
#                 .groups = "drop")
#     #merge the empirical accuracy back into the summaries
#     current.sim.summaries <- current.sim.summaries %>%
#       left_join(emp.accy, by = c("ModelName" = "ModelType", "Trait" = "fold"))
#     
#     #save it
#     write.table(current.sim.summaries, paste0("./Outputs/Summaries/Summary_", combo.tag, ".csv"), row.names = F)
#     
#     #write.csv(currrent.sim.blups, paste0("./Predictions/Preds_", combo.tag, ".csv"), row.names = F)
#     
#     gc()
#     
#    message(paste0("=== SUCCESS: COMPLETED JOB ", TASK_ID, " ==="))
