# MasterUtilityFunctions
# Sarah Lipps
# 2026-05-06
# These utility functions are for AIM1 - Optimizing use of secondary traits

######### Simulate Genomes and Traits #############
SimulateTraitsMarkers <- function(SP, POP, #GENOMES,
                                  rG12 = 0.7, rG13 = 0.1, rG23 = 0.0, 
                                  rE12 = 0.2, rE13 = 0.8, rE23 = 0.1,
                                  h2T1 = 0.2, h2T2 = 0.8, h2T3 = 0.3,
                                  nreps = 3, MISSING.PROP = 0) {
  ### CHECK THAT SP IS OF CLASS "MapPop" ###
  if(missing(SP)) {
    stop("Argument 'SP' is required and must be a 'SimParam' object generated with AlphaSimR \n")
  }
  if(!inherits(SP, "SimParam")) {
    stop("Argument 'SP' must be a 'SimParam' object generated using AlphaSimR \n")
  }
  
  print("Input objects have passed initial scan and validation! Proceeding with the simulation \n")
  
  #SP$traits <- list()
  if (SP$nTraits > 0) {
    SP$removeTrait(1:SP$nTraits)
  }
  SP$resetPed()
  SP$resetPed()
  
  ### GENOTYPE PARAMETERS ###
  # matrix for genetic correlation
  corA = matrix(c(
    1.0, rG12, rG13, #Trait 1: YIELD
    rG12, 1.0, rG23, #Trait 2: STANDARD 2o
    rG13, rG23, 1.0  #Trait 3: LATENT 2o
  ), nrow = 3, ncol = 3)
  
  # add the traits assuming 100 qtl/chr
  # qtl overlap is determined by rG
  SP$addTraitA(nQtlPerChr = 300, #maybe scale up to 300/chr
               mean = c(0, 0, 0),
               var = c(1, 1, 1),
               corA = corA)
  
  POP <- resetPop(POP, simParam = SP) # must come after addTraitA and before setPheno
  #SP$addSnpChip(nSnpPerChr = 500)
  #add the restrict sites...SnpChip doesn't overlap with QTL
  
  # create the population
  #pop = newPop(GENOMES)
  
  ### PHENOTYPE PARAMETERS ###
  corE = matrix(c(
    1.0, rE12, rE13,
    rE12, 1.0, rE23,
    rE13, rE23, 1.0
  ), nrow = 3, ncol = 3)
  
  # heritabilities for each trait
  my.h2s = c(h2T1, h2T2, h2T3)
  
  #Make repeated replications of the individuals and phenotypes
  nreps = nreps
  rep.list <- list()

  for(rep in 1:nreps){
    pop <- setPheno(POP, my.h2s, corE = corE, simParam = SP)
    # extract data values for trait
    tempdf <- data.frame(
      id = as.factor(pop@id),
      rep = as.factor(rep),
      Yield = pheno(pop)[, 1],
      Standard2o = pheno(pop)[, 2],
      Latent2o = pheno(pop)[, 3]
    )
    #store it
    rep.list[[rep]] <- tempdf
  }

  #make the df
  pheno.data <- do.call(rbind, rep.list)
  
  if(MISSING.PROP > 0){
    all.ids <- unique(pheno.data$id)
    ids2drop <- sample(all.ids, size = round(length(all.ids) * MISSING.PROP))
    pheno.data$Yield[pheno.data$id %in% ids2drop] <- NA
  }
  #results = list(pop, pheno.data)
  #return(results)
  
  
  message(" --> Trait Simulation Done <-- ")
  return(list(pop = POP, pheno.df = pheno.data))
}
 
######### Make Sparse Inverse Genomic Relationship Matrix #############
# Approach is Van Radden's as taught by Isik, Holland, & Malteca
# Additional code also from D. Tolhurst
MakeGRM <- function(SNP.MATRIX) {
  #adjusted to not require ASRgenomics
  #check that it is a matrix
  if(!is.matrix(SNP.MATRIX)) {
    stop("Argument 'SNP.Matrix' must be a matrix. You provided a ", class(SNP.Matrix))
  }
  #check that matrix is numeric
  if(!is.numeric(SNP.MATRIX)) {
    stop("Argument 'SNP.Matrix must be a numeric matrix.")
  }
  
  print("The input matrix has passed initial scan and validation! Proceeding with making your GRM!")
  
  M <- scale(SNP.MATRIX, scale = F)
  G <- M %*% t(M)
  G <- G/mean(diag(G))
  diag(G) <- diag(G) + 0.00001
  Ginv <- solve(G)
  
  I <- diag(nrow(G))
  G_blended <- 0.95 * G + 0.05 * I
  Ginv_stable <- solve(G_blended)
  # above avoids using ASR genomics
  message("The mean of the diagonal of you GENOMIC RELATION MATRIX is: ", mean(diag(G)))
  
  #adding to make sure it has rownames
  attr(Ginv_stable, "rowNames") <- rownames(G)
  attr(Ginv_stable, "INVERSE") <- TRUE
  return(list(G.inv = Ginv_stable, GRM = G))
}
######### Fit ASREML Model #############
# It is spaghetti format but it works
# Absolutely, under no circumstances DO NOT EDIT
m.asr <- function(GRM = NULL, PHENO.DAT, TRAITS.LIST, OPTION, USE.GRM = FALSE) {
  #check the option is valid
  valid.opts <- c("uni", "u", "biv", "b", "tri", "t", "help", "h")
  if(!(OPTION %in% valid.opts)) {
    warning(paste0("Option ", OPTION, " is invalid. Specify: uni/u, biv/b, tri/t, help/h."))
  }
  
  #help menu for options
  if(OPTION == "help") {
    cat("\n -- M.ASR OPTION GUIDE -- \n",
        "uni or u: Univariate for yield \n",
        "biv or b: Bivariate for yield, standard 2o or latent 2o \n",
        "tri or t: Trivariate for yield, standard 2o, and latent 2o \n",
        "help or h: help options \n")
    return(invisible(NULL))
  }
  
  if(USE.GRM & is.null(GRM)) {
    stop("You set USE.GRM = TRUE but did not provide a GRM matrix!")
  }
  
  #For storage and returning what is stored!
  mod.out.list <- list()
  pred.out.list <- list()
  
  #Run the ASREML model depending on the option you specify...
  if (OPTION == "uni" | OPTION == "u") {
    message("--- Running Univariate Model. GRM: ", USE.GRM, " ---")
    
    #loop through the traits provided
    for (trait in TRAITS.LIST) {
      fixed.part <- as.formula(paste0(trait, "~1+rep"))
      
      #determine the random part
      # if 2 then ~vm(id, GRM)
      # if 1 then ~id
      if(USE.GRM) {
        attr(GRM, "INVERSE") <- TRUE
        assign("GRM2", GRM, envir = .GlobalEnv)
        
        random.part <- as.formula(paste0("~vm(id,GRM2)"))
      } else {
        random.part <- as.formula("~id")
      }
      
      #run the model
      model <- tryCatch({
        asreml_args <- list(
          fixed = fixed.part,
          random = random.part,
          residual = as.formula(~units),
          data = PHENO.DAT,
          na.action = na.method(y = "include", x = "include"),
          workspace = '16gb',
          maxiter = 50,
          trace = FALSE
        )
        #update.count <- 0
        #max.updates <- 5
        temp.mod <- do.call(asreml, asreml_args)
        temp.mod
      }, error = function(e) {
        warning("ASReml FATAL ERROR: ", e$message)
        return(NULL)
      })
      #   while(!temp.mod$converge && update.count < max.updates) {
      #     message(" -> Model didn't converge. Updating...(Atempt", update.count + 1, ")")
      #     temp.mod <- update(temp.mod)
      #     update.count <- update.count + 1
      #   }
      #   temp.mod
      # }, error = function(e){
      #   warning("ASReml FATAL ERROR on trait: ", e$message)
      #   return(NULL)
      # })
      
      if(!is.null(model)) {
        message("...Getting the BLUPs! Hold tight!")
        mod.out.list[[trait]] <- model
        pred.out.list[[trait]] <- predict(model, classify = "id", only = "id",
                                          sed = TRUE, vcov = TRUE, pworkspace = '16gb')
      }
    } 
  } 
  
  if (OPTION == "biv" | OPTION == "b") {
    message("--- Running Bivariate Model. GRM: ", USE.GRM, " ---")
    
    for (trait.pr in TRAITS.LIST) {
      #Check that 2 traits are provided
      if (length(trait.pr) != 2) {
        warning(paste(" -> Skipping", trait.pr, "- not a pair of traits."))
        next
      }
      #t.name <- paste(trait.pr, collapse = "-")
      t1 <- trait.pr[1]
      t2 <- trait.pr[2]
      pair.name <- paste0(t1, "-", t2)
      
      fixed.part <- as.formula(paste0("cbind(",t1, ",", t2, ")~trait+trait:rep"))
      #determine the random part
      if(USE.GRM) {
        attr(GRM, "INVERSE") <- TRUE
        assign("GRM2", GRM, envir = .GlobalEnv)
        random.part <- as.formula("~us(trait):vm(id,GRM2)")
      } else {
        random.part <- as.formula("~us(trait):id")
      }
      residual <- as.formula("~units:us(trait)")
      #run the model
      model <- tryCatch({
        asreml_args <- list(
          fixed = fixed.part,
          random = random.part,
          residual = residual,
          data = PHENO.DAT,
          na.action = na.method(y = "include", x = "include"),
          workspace = '16gb',
          maxiter = 50,
          trace = F
        )
        #update.count <- 0
        #max.updates <- 5
        temp.mod <- do.call(asreml, asreml_args)
        temp.mod
      }, error = function(e) {
        warning("ASReml FATAL ERROR: ", e$message)
        return(NULL)
      })
      
      if (!is.null(model)) {
        message("...Getting the BLUPs! Hold tight!")
        mod.out.list[[pair.name]] <- model
        pred.out.list[[pair.name]] <- predict(model, classify = "trait:id", only = "trait:id",
                                           levels = list(trait=trait.pr[1]),
                                           sed = T, vcov = T, pworkspace = '16gb')
      }
    }
  } 
  if (OPTION == "tri" | OPTION == "t") {
    message("--- Running Trivariate Model. GRM: ", USE.GRM, " ---")
    
    for (trait.trio in TRAITS.LIST) {
      #check that 3 traits are provided
      if (length(trait.trio) != 3) {
        warning(paste(" -> Skipping", trait.trio, " - not a trio of traits."))
        next
      }
      t1 <- trait.trio[1]
      t2 <- trait.trio[2]
      t3 <- trait.trio[3]
      trio.name <- paste0(t1,"-", t2, "-", t3)
      fixed.part <- as.formula(paste0("cbind(",t1, ",", t2, ",", t3, ")~trait+trait:rep"))
      #determine the random part
      if(USE.GRM) {
        attr(GRM, "INVERSE") <- TRUE 
        assign("GRM2", GRM, envir = .GlobalEnv)
        random.part <- as.formula("~us(trait):vm(id,GRM2)")
      } else{
        random.part <- as.formula("~us(trait):id")
      }
      residual <- as.formula("~units:us(trait)")
      #run the model
      model <- tryCatch({
        asreml_args <- list(
          fixed = fixed.part,
          random = random.part,
          residual = residual,
          data = PHENO.DAT,
          na.action = na.method(y = "include", x = "include"),
          workspace = '16gb',
          maxiter = 50,
          trace = F
        )
        #update.count <- 0
        #max.updates <- 5
        temp.mod <- do.call(asreml, asreml_args)
        temp.mod
      }, error = function(e) {
        warning("ASReml FATAL ERROR: ", e$message)
        return(NULL)
      })
      
      if (!is.null(model)) {
        message("...Getting the BLUPs! Hold tight!")
        mod.out.list[[trio.name]] <- model
        pred.out.list[[trio.name]] <- predict(model, classify = "trait:id", only = "trait:id",
                                              levels = list(trait=trait.trio[1]),
                                              sed = T, vcov = T, pworkspace = '16gb')
      }
    } #end for loop
  } #end if statement
  
  return(list(models = mod.out.list, predictions = pred.out.list)) 
} #end function

######### Extract all of my Statistics #############
ExtractUni <- function(asreml.output, USE.GRM = FALSE, GRM = NULL, MODEL.NAME, SIM.NUM) {
  if(length(asreml.output$models) == 0) {
    return(list(
      summary.df = data.frame(
        SimulationNum = SIM.NUM, ModelName = MODEL.NAME, Trait = NA, Trait1 = NA,
        Trait2 = NA, Trait3 = NA, ModelType = ifelse(USE.GRM, "GBLUP", "StandardBLUP"),
        Vg1 = NA, Ve1 = NA, Vg2 = NA, Ve2 = NA, Vg3 = NA, Ve3 = NA,
        rg_12 = NA, re_12 = NA, rg_13 = NA, re_13 = NA,
        h2Standard = NA, h2Cullis = NA, Accuracy = NA,
        MeanPEV = NA, MeanAVSED = NA, Converged = FALSE, AIC = NA, BIC = NA
      ),
      blups = data.frame(
        id = character(), predicted.value = numeric(), PEV = numeric(),
        Gii = numeric(), r2 = numeric(), trait = character(),
        ModelType = character(), fold = character(), SimulationNum = numeric()
      )
    ))
  }
  #extract the named lists from the m.asr function output
  models.list <- asreml.output$models
  preds.list <- asreml.output$predictions
  
  #names of the traits that were run
  trait.names <- names(models.list)
  
  #empty storage containers (lists)
  summary.list <- list()
  updated.blups.list <- list()
  
  for (trait in trait.names) {
    
    #grab the specific model and prediction object for the list
    model <- models.list[[trait]]
    mod.summary <- summary(model)
    is.converged <- model$converge
    mod.aic <- mod.summary$aic
    mod.bic <- mod.summary$bic
    preds <- preds.list[[trait]]
    
    #Get the variance parameters (unscaled)
    vars <- model$sigma2 * model$vparameters
    
    if (USE.GRM) {
      Vg <- vars[grep("vm\\(id", names(vars))]
    } else {
      Vg <- vars[grep("^id$", names(vars))]
    }
    
    Ve <- vars[grep("units!R|residual", names(vars))]
    #if(is.na(Ve)) Ve <- model$sigma2
    if(length(Ve) == 0) Ve <- model$sigma2 
    
    #Calculate my metrics
    h2 <- Vg / (Vg + Ve) #plot level
    # add line mean here (Vg + Ve) / 3
    avsed <- preds$avsed[2]^2 / 2
    h2_Cullis <- 1 - (avsed / Vg)
    accy <- ifelse(!is.na(h2_Cullis) && h2_Cullis > 0, sqrt(h2_Cullis),  NA)
    
    #Looking at the individual predictions
    blups <- preds$pvals
    blups$PEV <- diag(preds$vcov)
    
    if (USE.GRM && !is.null(GRM)) {
      G.diag <- diag(GRM)
      blups$Gii <- G.diag[as.character(blups$id)]
      blups$r2 <- 1 - (blups$PEV / (blups$Gii * Vg))
    } else {
      blups$Gii <- 1
      blups$r2 <- 1 - (blups$PEV / Vg)
    }
    
    # Save the statistics as a 1 row summary
    summary.list[[trait]] <- data.frame(
      SimulationNum = SIM.NUM,
      ModelName = MODEL.NAME,
      Trait = trait,
      Trait1 = trait,
      Trait2 = NA,
      Trait3 = NA,
      ModelType = ifelse(USE.GRM, "GBLUP", "StandardBLUP"),
      Vg1 = Vg,
      Ve1 = Ve,
      Vg2 = NA,
      Ve2 = NA,
      Vg3 = NA,
      Ve3 = NA,
      rg_12 = NA,
      re_12 = NA,
      rg_13 = NA,
      re_13 = NA,
      h2Standard = h2,
      h2Cullis = h2_Cullis,
      Accuracy = accy,
      MeanPEV = mean(diag(preds$vcov)),
      MeanAVSED = avsed,
      Converged = is.converged,
      AIC = mod.aic,
      BIC = mod.bic
    )
    blups$trait <- trait #must include to match multivariate
    blups$ModelType = MODEL.NAME
    blups$fold <- trait
    blups$SimulationNum = SIM.NUM
    updated.blups.list[[trait]] <- blups
  }
  
  #bind all 1 row summaries into a single dataframe
  all.blups <- do.call(rbind, updated.blups.list)
  output.df <- do.call(rbind, summary.list)
  rownames(output.df) <- NULL
  rownames(all.blups) <- NULL
  #return everything
  return(list(
    summary.df = output.df,
    blups = all.blups
  ))
}
# Multivariate statistic extraction
ExtractMulti <- function(asreml.output, USE.GRM = FALSE, GRM = NULL, MODEL.NAME, SIM.NUM) {
  if(length(asreml.output$models) == 0) {
    return(list(
      summary.df = data.frame(
        SimulationNum = SIM.NUM, ModelName = MODEL.NAME, Trait = NA, Trait1 = NA,
        Trait2 = NA, Trait3 = NA, ModelType = ifelse(USE.GRM, "GBLUP", "StandardBLUP"),
        Vg1 = NA, Ve1 = NA, Vg2 = NA, Ve2 = NA, Vg3 = NA, Ve3 = NA,
        rg_12 = NA, re_12 = NA, rg_13 = NA, re_13 = NA,
        h2Standard = NA, h2Cullis = NA, Accuracy = NA,
        MeanPEV = NA, MeanAVSED = NA, Converged = FALSE, AIC = NA, BIC = NA
      ),
      blups =	data.frame(
        id = character(), predicted.value = numeric(), PEV = numeric(),
        Gii = numeric(), r2 = numeric(), trait = character(),
        ModelType = character(), fold = character(), SimulationNum = numeric()
      )
    ))
  }
  
  #get named list
  models.list <- asreml.output$models
  preds.list <- asreml.output$predictions
  #names of traits
  #trait.names <- names(models.list)
  #empty storage containers (lists)
  summary.list <- list()
  updated.blups.list <- list()
  
  for (trait in names(models.list)) {
    
    model <- models.list[[trait]]
    mod.summary <- summary(model)
    is.converged <- model$converge
    mod.aic <- mod.summary$aic
    mod.bic <- mod.summary$bic
    preds <- preds.list[[trait]]
    vars <- model$sigma2 * model$vparameters #double check this for the bivariate
    
    #split the names
    trait.elements <- unlist(strsplit(trait, "-"))
    t1 <- trait.elements[1]
    
    #tracking variables
    rg_12 <- NA; re_12 <- NA; t2 <- NA; Vg2 <- NA; Ve2 <- NA
    rg_13 <- NA; re_13 <- NA; t3 <- NA; Vg3 <- NA; Ve3 <- NA
    
    # Extract Vg and covg for trait 1 
    if (USE.GRM) {
      Vg1 <- vars[grep(paste0("vm.*trait_", t1, ":", t1), names(vars))] 
      #Vg1 <- vars[grep(paste0("vm.*!trait_", t1, ":", t1), names(vars))]
    } else {
      Vg1 <- vars[grep(paste0("trait:id.*trait_",t1,":",t1), names(vars))] 
      #Vg1 <- vars[grep(paste0("trait:id!trait_", t1, ":", t1), names(vars))]
    }
    #get VE for trait 1
    Ve1 <- vars[grep(paste0("units:trait!trait_", t1, ":", t1), names(vars))]
    if(length(Ve1) == 0) Ve1 <- model$sigma2
    
    
    #if bivariate (2 traits )
    if(length(trait.elements) >= 2) {
      t2 <- trait.elements[2]
      if (USE.GRM) {
        Vg2 <- vars[grep(paste0("vm.*!trait_", t2, ":", t2), names(vars))]
        COVg12 <- vars[grep(paste0("vm.*!trait_(", t1, ":", t2, "|", t2, ":", t1, ")"), names(vars))]
      } else {
        Vg2 <- vars[grep(paste0("trait:id!trait_", t2, ":", t2), names(vars))]
        COVg12 <- vars[grep(paste0("trait:id!trait_(", t1, ":", t2, "|", t2, ":", t1, ")"), names(vars))]
      }
      if(length(Vg1) > 0 && length(Vg2) > 0 && length(COVg12) > 0) {
        rg_12 <- COVg12 / sqrt(Vg1 * Vg2)
      }
      #residual
      Ve2 <- vars[grep(paste0("units:trait!trait_", t2, ":", t2), names(vars))]
      COVe12 <- vars[grep(paste0("units:trait!trait_(", t1, ":", t2, "|", t2, ":", t1, ")"), names(vars))]
      if(length(Ve1) > 0 && length(Ve2) > 0 && length(COVe12) > 0) {
        re_12 <- COVe12 / sqrt(Ve1 * Ve2)
      }
    }
    
    #if traivariate (3 traits)
    if(length(trait.elements) >= 3) {
      t3 <- trait.elements[3]
      if(USE.GRM) {
        Vg3 <- vars[grep(paste0("vm.*!trait_", t3, ":", t3), names(vars))]
        COVg13 <- vars[grep(paste0("vm.*!trait_(", t1, ":", t3, "|", t3, ":", t1, ")"), names(vars))]
      } else {
        Vg3 <- vars[grep(paste0("trait:id!trait_", t3, ":", t3), names(vars))]
        COVg13 <- vars[grep(paste0("trait:id!trait_(", t1, ":", t3, "|", t3, ":", t1, ")"), names(vars))]
      }
      if(length(Vg1) > 0 && length(Vg3) > 0 && length(COVg13) > 0) {
        rg_13 <- COVg13 / sqrt(Vg1 * Vg3)
      }
      #residual
      Ve3 <- vars[grep(paste0("units:trait!trait_", t3, ":", t3), names(vars))]
      COVe13 <- vars[grep(paste0("units:trait!trait_(", t1, ":", t3, "|", t3, ":", t1, ")"), names(vars))]
      if(length(Ve1) > 0 && length(Ve3) > 0 && length(COVe13) > 0) {
        re_13 <- COVe13 / sqrt(Ve1 * Ve3)
      }
    } 
    
    #calculate metrics for trait 1
    h2 <- Vg1 / (Vg1 + Ve1)
    avsed <- preds$avsed[2]^2/2
    h2_Cullis <- 1 - (avsed / Vg1) #this is actually (generalized) reliability b/c on an estimator
    accy <- ifelse(!is.na(h2_Cullis) && h2_Cullis > 0, sqrt(h2_Cullis), NA)
    
    #get the individual predictions and reliability
    blups <- preds$pvals
    blups$PEV <- diag(preds$vcov)
    
    if(USE.GRM && !is.null(GRM)) {
      G.diag <- diag(GRM)
      blups$Gii <- G.diag[as.character(blups$id)]
      blups$r2 <- 1 - (blups$PEV / (blups$Gii * Vg1))
    } else {
      blups$Gii <- 1
      blups$r2 <- 1 - (blups$PEV / Vg1)
    }
    
    #save the 1-row summary
    summary.list[[trait]] <- data.frame(
      SimulationNum = SIM.NUM,
      ModelName = MODEL.NAME,
      Trait = trait,
      Trait1 = t1,
      Trait2 = t2,
      Trait3 = t3,
      ModelType = ifelse(USE.GRM, "GBLUP", "StandardBLUP"),
      Vg1 = Vg1,
      Ve1 = Ve1,
      Vg2 = Vg2,
      Ve2 = Ve2,
      Vg3 = Vg3,
      Ve3 = Ve3,
      rg_12 = rg_12,
      re_12 = re_12,
      rg_13 = rg_13,
      re_13 = re_13,
      h2Standard = h2,
      h2Cullis = h2_Cullis,
      Accuracy = accy,
      MeanPEV = mean(diag(preds$vcov), na.rm = T),
      MeanAVSED = avsed,
      Converged = is.converged,
      AIC = mod.aic,
      BIC = mod.bic
    )
    
    blups$ModelType = MODEL.NAME
    blups$fold <- trait
    blups$SimulationNum <- SIM.NUM
    updated.blups.list[[trait]] <- blups
  }
  all.blups <- do.call(rbind, updated.blups.list)
  output.df <- do.call(rbind, summary.list)
  rownames(output.df) <- NULL
  rownames(all.blups) <- NULL
  return(list(
    summary.df = output.df,
    blups = all.blups
  ))
}
