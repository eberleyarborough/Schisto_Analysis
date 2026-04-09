####
# INLA Helper Functions
# Created 2025 - adapted from HJC code
####

### Data subset function
my_dt_subset_function <- function(DT = NULL, Covar = NULL, Resp = NULL){
  dt.sub <- DT[complete.cases(DT[,colnames(DT) %in% Covar]),] 
  dt.sub <- dt.sub[complete.cases(dt.sub[,colnames(dt.sub) %in% Resp]),]  
  dt.sub
}
  
  
### Location format function
my_location_function <- function(DT = NULL, LOCUTM = FALSE){
  if(isFALSE(LOCUTM)){
    require(sp)
    
    c.utm <- as.data.frame(cbind(DT$x.east, DT$y.south))
    names(c.utm) <-  c("X","Y")
    
    coordinates(c.utm) <- c("X","Y")
    proj4string(c.utm) <- CRS("+proj=longlat +datum=WGS84")
    
    c.utm <- spTransform(c.utm, CRS("+proj=utm +zone=36 ellps=WGS84"))
    
    # Extracting the coordinates themselves
    LOC <- coordinates(c.utm)  
    
  } else {
    LOC <- cbind(DT$x.east, DT$y.south)
  }
  
  return(LOC)
}


####
# Mesh function
####
my_mesh_function <- function(DT = NULL,          # Data frame with observations
                             Covar = NULL,       # Covariate names (for subsetting)
                             Resp = NULL,        # Response variable name (for subsetting)
                             LOCUTM = FALSE,     # Are coordinates in UTM or do they need conversion? (T/F)
                             toUTM = NULL,       # Currently unused parameter, not sure what Henri did with this
                             BNDRY = NULL,       # Boundary object defining the study area limits (here is `kruger_shape`)
                             Max_edge = NULL,    # Maximum triangle edge length [inner, outer], set previously
                             Cutoff = NULL,      # Minimum distance between edge nodes, set previously
                             Offset = NULL,      # Extension beyond boundary [inner, outer], set previously
                             DTsubset = FALSE){  # Should the data be subset to complete cases? (optional)
  
  # ===== STEP 1: OPTIONAL DATA SUBSETTING =====
  # If DTsubset=TRUE, remove rows with missing values in covariates or response (presumably we've already done this)
  # Actual subsetting is taken care of with the `my_dt_subset_function` code
  # This ensures the mesh is built only for locations with complete data
  if(isTRUE(DTsubset)){
    # Removes NAs so mesh nodes are only placed where we have usable data
    DT <- my_dt_subset_function(DT, Covar = Covar, Resp = Resp)
  } 
  
  # ===== STEP 2: EXTRACT LOCATION COORDINATES =====
  # Get X,Y coordinates for all observations
  # These are the "knots" - actual data points where observations were made
  # The mesh will be constrained to cover these locations plus surrounding area within KNP
  LOC <- my_location_function(DT = DT, LOCUTM = LOCUTM)
  
  # Debug output to verify coordinates are in the expected range
  # For KNP UTM: X should be ~280k - 400k, Y should be ~7.1M - 7.5M
  cat("LOC in mesh function - X range:", range(LOC[,1]), "\n")
  cat("LOC in mesh function - Y range:", range(LOC[,2]), "\n")
  
  # ===== STEP 3: VERIFY BOUNDARY COORDINATES =====
  # Check that the boundary coordinates match the data coordinate system to make sure everything is projecting correctly
  # Mismatched coordinate systems will create an invalid mesh
  if(!is.null(BNDRY)) {
    cat("BNDRY$loc - X range:", range(BNDRY$loc[,1]), "\n")
    cat("BNDRY$loc - Y range:", range(BNDRY$loc[,2]), "\n")
  }
  
  # ===== STEP 4: CREATE THE MESH =====
  # This is the core function that builds the triangular mesh
  # NOTE: fmesher package replaced the old INLA mesh functions
  # fm_mesh_2d_inla is the modern equivalent of inla.mesh.2d
  Mesh = fmesher::fm_mesh_2d_inla(
    loc = LOC,                               # Observation coordinates (n x 2 matrix)
    boundary = fmesher::fm_as_segm(BNDRY),   # Convert boundary to fmesher format
    max.edge = Max_edge,                     # Triangle size control [inner, outer]
    cutoff = Cutoff,                         # Minimum node spacing
    offset = Offset)                         # Boundary extensions [inner, outer]
  
  # ===== STEP 5: VERIFY MESH COORDINATES =====
  # Check that mesh was created in the correct coordinate system
  # Mesh coordinates should match LOC and BNDRY coordinates
  cat("Created mesh - X range:", range(Mesh$loc[,1]), "\n")
  cat("Created mesh - Y range:", range(Mesh$loc[,2]), "\n")
  
  # ===== STEP 6: RETURN THE MESH =====
  # The mesh object contains:
  # $loc, the coordinates of all mesh vertices (nodes)
  # $graphh, the triangulation structure (which nodes connect to which)
  # $n, the total number of mesh nodes
  # This mesh will be used to create the SPDE model and projection matrix (A)
  Mesh
}


####
# spde function - here is where we set SPDE priors and controls
###

my_spde_function <- function(Mesh, spdePriorFlat = TRUE,
                             spdepriorrange = NULL,
                             spdepriorsigma = NULL){
  
  # ===== OPTION 1: FLAT (VAGUE) PRIORS =====
  # Uses default INLA priors - less informative, lets data drive estimates
  # Good for when you have no prior knowledge about spatial patterns that might be in the system
  # Drawback can lead to overfitting or unrealistic spatial range estimates  
  if(isTRUE(spdePriorFlat)){
    buff.Hosts.spde = inla.spde2.matern(mesh = Mesh)
    # Creates Matérn SPDE model with default priors
    # Matérn covariance is the standard choice for spatial modeling
    # Assumes spatial correlation decays smoothly with distance
  }
  
  # ===== OPTION 2: PENALIZED COMPLEXITY (PC) PRIORS =====
  # Uses informative priors based on prior knowledge of the study system
  # PC priors penalize complexity - in stats, it's good practice to prefer simpler models unless data strongly suggests otherwise
  # Good for when you have scientific knowledge about the spatial scale of whatever processes you're trying to model
  # This bit of code says what to do if the spdePriorFlat is FALSE - use the priors we set earlier in our INLA list setup!
  else {
    buff.Hosts.spde = inla.spde2.pcmatern(
      mesh = Mesh, 
      prior.range = spdepriorrange,     # Prior for spatial range parameter (spatial correlation decay distance)
      prior.sigma = spdepriorsigma)     # Prior for field variance parameter (marginal std. dev of spatial field, aka size of variability)
  }
  
  return(buff.Hosts.spde)
}


#
# FUNCTION: my_modelmatrix
#
# used in Stack generation below
my_modelmatrix_function <- function(DT = NULL, Covar = NULL, CatVarSub = NULL){
  
  # Convert any factors or integers to characters (so dummy vars are created in the matrix) 
  # Note: this creates ONE, BINARY covariate per each character level - e.g., SexMale, SexFemale for `Sex`, which has 2 levels
  # This is because in Bayesian stats, it's best practice to include a column for each level, with EACH representing that level's deviation from the global mean (i.e., the intercept)
  # E.g., `SexMale` is the deviation from the mean associated with males, and `SexFemale` is the deviation from the overall mean associated with females
  DT[,Covar] %>% dplyr::mutate(across(where(is.factor), as.character),
                               across(where(is.integer), as.character)) -> DT.modmatrix
  
  X0 <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(Covar, collapse = " + "))),
                     data = DT.modmatrix)
  
  # Remove columns that are literally "(Intercept)" if they exist
  X <- as.data.frame(X0[, !colnames(X0) %in% "(Intercept)"])
  X
}

# Testing the output to make sure no covariates are dropped
#X_test <- my_modelmatrix_function(DT = positive_only_schisto_df, Covar = covar_mod2)
#colnames(X_test)


####
# Stack function
####

my_stack_function <- function(
    DT = NULL, 
    Mod.Matrix = NULL, 
    Covar = NULL, 
    Resp = NULL, 
    GroupVar = NULL, 
    spdePriorFlat = FALSE, 
    spdepriorrange = NULL, 
    spdepriorsigma = NULL, 
    LOCUTM = NULL, 
    Mesh.list = NULL, 
    BNDRY = NULL,  
    Max_edge = Max_edge_set, 
    Predict = FALSE,         # Where you update when you want the model to generate the prediction outputs
    include_section = NULL,  # Where you'd want to update if your random effect (in this case, section) changes
    include_spatial = FALSE){ # Default is set to false so spde elements don't get included in models 1 and 2
  
  # ===== STEP 1: GET LOCATION COORDINATES =====
  # Extract spatial coordinates (X, Y) from data
  # If LOCUTM=TRUE, transforms lat/long to UTM; if FALSE, uses existing coordinates
  LOC <- my_location_function(DT = DT, LOCUTM = LOCUTM)
  
  # ===== OPTIONAL: CREATE ID VARIABLES FOR RANDOM EFFECTS =====
  # See Henry's code `_hjc` for an example
  
  # ===== STEP 3: PREPARE THE COVARIATE MATRIX =====
  # Convert the model matrix "X" (which contains covariates like Sex, BCS, age, etc.) from a matrix into a list of individual named columns
  # This allows INLA to reference covariates by name (e.g., "SexFemale")
  # The model matrix was created outside this function in the `my_modelmatrix_function`
  if(!is.null(Mod.Matrix)) {
    X_list <- as.list(as.data.frame(Mod.Matrix))
    names(X_list) <- colnames(Mod.Matrix)
  } else {
    X_list <- list()
  }
  
  # ===== STEP 4: BUILD EFFECTS LIST FOR ALL MODELS =====
  # Here we put all of the fixed effects into one list
  # If a variable is in the formula, it must live in effects, not data
  fixed_effects <- c(
    list(Intercept = rep(1, nrow(DT))),
    X_list
  )
  
  # Start with fixed effects as FIRST element
  # Effects list = what lives where
  effects_list <- list(fixed_effects)
  
  # ===== STEP 6: INITIALIZE A LISTS =====
  # A list = how the things in the effects list are projected onto the observations
  A_list <- list(1)  # ONE projection for ALL fixed effects
  
  # ===== STEP 5: (OPTIONAL) ADD HERD-LEVEL RANDOM EFFECT =====
  # This code is set up to include park section as a random effect
  if(include_section) {
    NSection <- as.numeric(as.factor(DT$Section))
    effects_list <- c(effects_list, list(list(section.num = NSection)))  # SECOND element
    A_list <- c(A_list, list(1))  # SECOND identity projection for section
  }
  
  # ===== STEP 7: ADD SPATIAL COMPONENT (ONLY IF REQUESTED)
  # This is where spatial random effects are added via SPDE
  if(include_spatial) {
    
    # Get location coordinates again (redundant, but ensures correct coordinates)
    LOC <- my_location_function(DT = DT, LOCUTM = LOCUTM)
    Mesh = Mesh.list
    
    # ===== CREATE PROJECTION MATRIX A =====
    # buff.HostsA is the KEY matrix that projects observation locations onto the mesh (not to be confused with the A matrix earlier)
    # It's a sparse matrix: [n_observations × n_mesh_nodes]
    # Each row corresponds to one observation; columns are mesh vertices
    # Values indicate how much each mesh node contributes to that observation
    if(!is.null(GroupVar)){
      # If there's a grouping variable (e.g., year), create grouped spatial effects (we don't have this currently)
      buff.HostsA <- inla.spde.make.A(Mesh, 
                                      loc = LOC, 
                                      group = as.numeric(as.factor(DT[,GroupVar])))
    } else {
      # Standard case: just project locations onto mesh
      buff.HostsA <- inla.spde.make.A(Mesh, loc = LOC)
    }
    
    cat("buff.HostsA dimensions:", dim(buff.HostsA), "\n")
    
    # ===== CREATE SPDE MODEL =====
    # This defines the spatial correlation structure
    # Pulls in our previously specified priors for range (how far spatial  autocorrelation extends) and sigma (field variance)
    buff.Hosts.spde <- my_spde_function(Mesh,
                                        spdePriorFlat = spdePriorFlat,
                                        spdepriorrange = spdepriorrange,
                                        spdepriorsigma = spdepriorsigma)
    
    # ===== CREATE SPATIAL FIELD INDEX =====
    # s.index assigns an index number to each mesh node
    # This creates the "spatial.field" variable that can be used in formulas
    # Extract only the FIRST component (see below)
    s.index <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = buff.Hosts.spde$n.spde)
    
    # Extract only the spatial field indices (not group or replicate components) and add to effects list as another element
    effects_list <- c(effects_list, list(spatial.field = s.index$spatial.field))
    
    # Add the spatial projection matrix `buff.HostsA` to `A_list`
    # This is the actual spatial component of the projection matrix now that is being COMBINED with the identity projection list we made earlier
    A_list <- c(A_list, list(buff.HostsA))
  }
  
  # ===== STEP 9: DEBUG OUTPUT =====
  cat("\n=== INSIDE STACK FUNCTION - PRE-STACK ===\n")
  cat("nrow(DT):", nrow(DT), "\n")
  cat("length(y):", length(DT[,Resp]), "\n")
  #cat("\nA_list length:", length(A_list), "\n")
  cat("effects_list length:", length(effects_list), "\n")
  if(length(effects_list) > 0) {
    cat("effects_list names:", paste(names(effects_list), collapse=", "), "\n")
  }
  cat("=========================================\n\n")
  
  # ===== STEP 10: BUILD THE INLA STACK =====
  # The INLA stack packages together 3 things:
  # The data, or the response variable (y)
  # The A list, or the list of projection matrices (how to map effects to observations)
  # The effects, or a list of all predictor variables and indices
  # INLA uses this "stack" structure to understand the model's data organization
  
  stack_data <- list(y = DT[,Resp])  # the stack data contains ONLY the response variable, `y`
  
  stack_model <- inla.stack(
    data = stack_data,                # Only `y`
    A = A_list,                       # All projection matrices
    effects = effects_list,           # All identity projections (if included)
    tag = "model_stack"               # Name tag for this specific stack (as opposed to "predicted stack")
  )
  
  cat("\n=== DEBUGGING STACK INTERNALS ===\n")
  cat("Checking effects_list before stack creation:\n")
  
  cat("\n=== IMMEDIATELY AFTER STACK CREATION ===\n")
  test_data <- inla.stack.data(stack_model)
  cat("Test extraction - y length:", length(test_data$y), "\n")
  cat("Test extraction - Intercept length:", length(test_data$Intercept), "\n")
  cat("========================================\n\n")
  
  if (!isTRUE(Predict)) {
    return(stack_model)
  }
  
  # ===== STEP 11: CREATE THE PREDICTION STACK (OPTIONAL) =====
  # If Predict=TRUE, create a second stack for making predictions
  # This has the same structure but with y=NA (unknown values to predict)
  if(isTRUE(Predict)){
    
    # Build prediction effects list conditionally (same structure as model effects)
    pred_fixed_effects <- c(
      list(Intercept = rep(1, nrow(DT))),
      as.list(as.data.frame(Mod.Matrix))
    )
    
    pred_effects_list <- list(pred_fixed_effects)
    pred_A_list <- list(1)
    
    # Add section (or other random effect) if it was included in the model
    if(include_section) {
      pred_effects_list <- c(pred_effects_list, list(section.num = NSection))
      pred_A_list <- c(pred_A_list, list(1))  # Use list(1), not just 1 to give it it's own identity projection
    }
    
    # Add spatial component if it was included in the model
    if(include_spatial) {
      pred_effects_list <- c(pred_effects_list, list(spatial.field = s.index$spatial.field))
      pred_A_list <- c(pred_A_list, list(buff.HostsA))
    }
    
    # Create prediction stack with NA response values
    stack_predict <- inla.stack(
      data = list(y = NA),             # Unknown values to predict
      A = pred_A_list,
      effects = pred_effects_list,
      tag = "stack_predict"            # Different tag to distinguish from the model stack
    )
    
    # Combine model and prediction stacks
    # INLA will fit the model using `stack_model` data
    # Then use `stack_predict` structure to generate predictions
    return(inla.stack(stack_model, stack_predict))
  } 
}  


#### END CUSTOM FUNCTIONS ####
