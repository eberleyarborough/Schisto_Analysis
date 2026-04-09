########
#' This script is a variation of the `schisto_inlabru_v1.Rmd` document,
#' which seeks to understand drivers of Schistosoma. spp. infecion in African
#' buffalo within KNP. 

#' Preliminary work showed that perhaps there are very different spatial patterns
#' and/or drivers in the north vs south of KNP, with the most recent (3.9.26)
#' prevalence model finding a focal hotspot just south of Skukuza. It is possible
#' that this hotspot is masking any apparent patterns in the north of the park.
#' Therefore, we decided to split KNP into "North" and "South" regions, to
#' generate a flexible model to capture differences in covariate effects
#' between those regions of the park.

#' Previous models have illustrated that for prevalence, a Bernoulli GAM + SRF
#' might best model schisto in buffalo within KNP, so we will build off of that
#' model here. We can compare this model with flexible region components to our
#' previous model where buffalo across the entire park share 1 SRF.

#' We'll model with the following equation, which is the same as model 2 in the 
#' `schisto_inlabru_v1.Rmd` file, but our SRF here is split into two, independent 
#' SRFs by North and South:

#'    result01_ir ~ Bernoulli(Pi_ir)
#'    [result01] = Pi_ir
#'    var[result01] = Pi_ir * (1 - Pir)
#'    
#'    logit_(Pi_ir) = Intercept + Covariates + u_r(i), 
#'    where u_r(i) are the spatially correlated random effects that are 
#'    allowed to differ by region (observation i at region r)

#' We can write this logit link function as: 

#'          exp(Intercept + Covariates_i + u_r(i))  
#' Pi_i = ------------------------------------------
#'        1 + exp(Intercept + Covariates_i+ u_r(i))  

#' The SRFs for the two regions will be called SRF.s and SRF.n, and they will be
#' defined on the same mesh and use the same priors, but they will be estimated
#' separately for both south and north regions. You can choose different priors
#' and meshes for the two regions if you want to explore that.

########

# Section 1: Set up -----

#' We need the following packages:
library(lattice)
library(ggplot2)
library(sf)
library(gstat)
library(INLA)
library(inlabru)
library(DHARMa)
library(dplyr)
library(glmmTMB)
library(rnaturalearth)
library(scales)
library(cowplot)
library(tidyr)
source("/Users/eberleyarborough/Documents/Coding Help/INLA tutorials/highland_stats_INLA_course/AllRCodeInlabru/HighstatLibV14.R") #' Highland Statistics support file

#' Set the working directory (should be the project directory).
setwd("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/Schisto/Schisto_Analysis")



# Section 2: Data Importing & Exploration -----


#* Subsection 2.1: Importing the data -----

#' Step 1: Import the disease data
#' This file contains the schisto results, as well as the buffalo demographics,
#' capture locations, and elevation at those capture locations:
schisto_base<- read.csv("schisto_with_cap_demo_and_elev.csv")

#' converting to a spatial object
schisto.sf <- st_as_sf(schisto_base,
                       coords = c("x.east", "y.south"),
                       crs    = 32736)

#' Setting the CRS, converting to KM (it is best to work with small numbers in inlabru)
crs.Target <- fm_crs("EPSG:32736")  # Target for KNP in UTM 36S
CRS.km     <- fm_crs_set_lengthunit(crs.Target, "km") # Converting to km instead of m

#' We now transform the spatial object to the kilometer-based coordinate 
#' system and extract the coordinates as numeric variables. These will 
#' be used later when building spatial models.
schisto2.sf <- st_transform(schisto.sf, CRS.km)
Coords <- st_coordinates(schisto2.sf)
schisto2.sf$Xkm <- Coords[, 1]
schisto2.sf$Ykm <- Coords[, 2]

#' Dropping unnecessary columns from the schisto data frame
schisto2.sf <- schisto2.sf %>%
  select("animal.ID", "age.years", "age.class", "Sex", "BCS.Average", "Result", "Xkm", "Ykm", "geometry")

#' Step 2: Import the environmental covariates
#' Read in the .csv file with the environmental covariates. NOTE: This file
#' contains data for all capture locations, not just those that are schisto positive.
env_data <- read.csv("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/CSVs/cap_and_env_data.csv")
env_data <- env_data %>%
  select("animal.ID", "capture.date", "capture.day", "capture.group.by.day", "approx.herd.size", "weighted_dist_to_fire", "topo_wetness_index", "clay_percent_30cm", "jan_apr_23_evi", "jan_apr_23_msi")

#' Join to schisto data
schisto.complete.sf <- left_join(schisto2.sf, env_data, by = "animal.ID")

#' Step 3: Transform Data

#' It's good to convert the character "positive/negative" results to 
#' something binary and numeric (0/1, where 1 is positive)
schisto.complete.sf <- schisto.complete.sf %>%
  mutate(result01 = ifelse(Result == "positive", 1, 0))

#' And convert Sex, age class, and capture group by day to factors
schisto.complete.sf$fage.class <- factor(schisto.complete.sf$age.class)
schisto.complete.sf$fSex <- factor(schisto.complete.sf$Sex)
schisto.complete.sf$fcapgroup <- factor(schisto.complete.sf$capture.group.by.day)

#' We'll also log1p transform age in case we want to model with it later on. We
#' also use the `log1p()` function here instead of just `log()` because we would
#' end up zero-inflating our dataset by taking the ln of 1 ( = 0).
schisto.complete.sf <- schisto.complete.sf %>%
  mutate(log.age = log1p(age.years))

#' Step 4: Import the KNP boundary and match CRS
#' To provide geographic context, we load a map of Kruger that will be 
#' used to visualize the study area and the sampling locations.
kruger_shape <- st_read("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/SANParks KNP Shapefiles/KNP_polygon/KNP_polygon.shp")

#' Transform to UTM Zone 36S in km
kruger.km <- st_transform(kruger_shape, crs = CRS.km)

#' Finally, we plot the KNP boundary together with the spatial distribution of 
#' the capture locations to confirm that the data are correctly georeferenced. 
ggplot() +
  geom_sf(data = kruger.km, fill = "transparent") +
  geom_sf(data = schisto.complete.sf) +
  theme_minimal()

#' Confirm that the data and the map share the same CRS: 
st_crs(kruger.km) == st_crs(schisto.complete.sf)
#' output should be "TRUE"




#* Subsection 2.2: Exploring the data -----

#' Step 1: Determine if there are any missing values
MyVar <- c("animal.ID", "age.years", "age.class", "log.age", "Sex", "BCS.Average", "result01",
           "capture.date", "capture.day", "capture.group.by.day", "approx.herd.size",
           "weighted_dist_to_fire", "topo_wetness_index", "clay_percent_30cm",
           "jan_apr_23_evi", "jan_apr_23_msi")
colSums(is.na(schisto.complete.sf[,MyVar]))
#' Looks like missing values in herd size, EVI, and MSI (11 in each).

#' We'll drop NA rows
#' Now we'll remove any NA columns
schisto.complete.sf <- schisto.complete.sf %>%
  drop_na()
#' Left with 557 observations

#' For some analyses, it is convenient to work with a standard data frame 
#' rather than an sf object. We therefore drop the geometry column and 
#' store the result separately.
schisto.complete.df <- st_drop_geometry(schisto.complete.sf)

#' It's also a good idea to standardize continuous covariates - this can only be done after NA is dropped.
#' The `MyStd()` function is a custom function in the Highland Statistics 
#' script.
schisto.complete.sf$weighted_dist_to_fire.std <- MyStd(schisto.complete.sf$weighted_dist_to_fire)
schisto.complete.sf$topo_wetness_index.std <- MyStd(schisto.complete.sf$topo_wetness_index)
schisto.complete.sf$clay_percent_30cm.std <- MyStd(schisto.complete.sf$clay_percent_30cm)
schisto.complete.sf$jan_apr_23_evi.std <- MyStd(schisto.complete.sf$jan_apr_23_evi)
schisto.complete.sf$jan_apr_23_msi.std <- MyStd(schisto.complete.sf$jan_apr_23_msi)

#' Refer to `schisto_inlabru_v1.Rmd` for a more complete data exploration.

#' Step 2: Assess for collinearity
#' Make a pairplot - the `Mypairs()` function comes from the Highland Statistics
#' custom script.
MyVar <- c("Xkm", "Ykm", "log.age", "BCS.Average", "clay_percent_30cm",
           "jan_apr_23_evi", "jan_apr_23_msi", "topo_wetness_index",
           "weighted_dist_to_fire")
Mypairs(schisto.complete.df[,MyVar])
#' There might be some collinearity between Xkm and/or Ykm and percentage of clay
#' in the soil, as well as EVI and MSI, which makes sense, as these are spatial
#' covariates.

#' We only should be nervous if a Pearson correlation is larger than 0.8; however,
#' that doesn't mean we should ignore moderately high (~0.6) values, like
#' those for X/Ykm and clay percentage.

#' Step 3: Put the data into "long format" for easier plotting.
GatherTheseCols <- c("jan_apr_23_msi", "jan_apr_23_evi", "topo_wetness_index",
                     "weighted_dist_to_fire", "age.years", "log.age")
schist_long <- gather(data = schisto.complete.df,
                      key = "ID",
                      value = "AllX",
                      all_of(GatherTheseCols),
                      factor_key = TRUE)
dim(schist_long)

#' Step 3: Continuous covariates vs. disease data
#' Plot each continuous covariate against the presence/absence data. Note:
#' the disease data MUST be in numeric 0/1 form to see if relationships
#' are non-linear with this plot.
ggplot(data = schist_long, aes(x=AllX, y = result01)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Covariate", y = "Disease Status") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 5)) +
  geom_smooth(method = "gam",
              method.args = list(family = "binomial"),
              se = FALSE) +
  facet_wrap(.~ID, scale = 'free_x', ncol = 2)
#' Only non-linear pattern appears to be weighted distance to fire, so 
#' we might need a smoother for this covariate.


#* Subsection 2.3: Splitting the park -----

#' Here, we need to separate capture locations that occurred north and south
#' of Satara rest camp into "North" and "South" capture locations, respectively.
#' We'll then generate separate SRFs by region to include in a model. We happen
#' to have a shapefile from SANParks of all KNP camps, so we'll try to pull the 
#' Satara coordinates from that.

#' Read in the SANParks rest camp data
camps <- st_read("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/KNP Environmental Data/SANParks_data_15012026/camps_all.shp")

head(camps)
#' This shapefile contains the name, description, borehole, operational status, 
#' colloquial synonym, access type, type of camp, lat, long, status, and point location
#' in UTM for every camp in KNP.

#' Filter for Satara-specific information
satara <- camps %>%
  filter(NAME == "Satara")

#' Ensure that satara is in the right CRS and in KM
satara <- st_transform(satara, CRS.km)

#' Extract coordinates
satara.coords <- st_coordinates(satara)
satara.northing <- satara.coords[, "Y"]
satara.easting <- satara.coords[, "X"]

#' And group buffalo capture locations by North/South relationship to Satara
schisto.complete.sf <- schisto.complete.sf %>%
  mutate(region = ifelse(Ykm > satara.northing, "North", "South"))

#' See how many we have per region
table(schisto.complete.sf$region)
#' 368 in the north, 189 in the south

#' Out of curiosity, look at positives/negatives by north/south
table(schisto.complete.sf$region, schisto.complete.sf$result01)
#' We don't have complete separation, which is good. Overall, there is a higher
#' proportion of positives in the south compared to the north.


# Section 3: Building the Mesh -----

#' Step 1: Making the mesh -----

#' These settings are copy and paste from the `schisto_inlabru_v1.Rmd` file.
Loc <- cbind(schisto.complete.sf$Xkm, schisto.complete.sf$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
#' Sites are anywhere from 0-350km apart, so anything small scale is probably < 50 km.

#' The main tools to control the mesh are:
#' - max.edge: The largest allowed triangle length. The lower the value 
#'             for max.edge the higher the resolution. This should be "small-scale"
#'             based off of your sampling locations (1/5 the resolution).
#' - cutoff:   Sites with a distance value smaller than the cutoff 
#'             are modeled by a single vertex in the mesh.
#'             Set the cutoff to 1/5 of the maxedge of the inner area.
#' - boundary: Useful for islands and fjords, but you can make it hug your sampling
#'             locations tighter by making it nonconvex.
#  - offset:   Determines how far the inner and outer boundaries extend.
#' - Inner area / outer areas:
#'             Is to avoid potential numerical estimation problems for 
#'             the points close to the boundary. General recommendation is: 1:5.
#'             This explains the:   max.edge = c(1, 5) * MaxEdge in the 
#'             code below.

#' Now we should define boundaries for the spatial mesh.
#' First, define the resolution of the mesh based off of our distances between sites.
#' This can be adjusted to get at finer-scale (local) spatial patterns. 
MeshResolution <- 15               #' in km.
MaxEdge  <- MeshResolution / 5    #' good rule of thumb is to have max edge set to 1/5 the resolution

#' #' fm_nonconvex_hull makes a non-convex area around the sampling
#' locations - the inner boundary. You can control how close the blue line is 
#' to the sites with the 'convex' argument. 
inner.hull <- fm_nonconvex_hull(Loc, convex = -0.05, crs = st_crs(schisto.complete.sf))
#' Results
plot(inner.hull)
points(Loc)

#' Convert the inner, nonconvex hull to an sf object
inner.hull_sf <- st_as_sf(inner.hull)

#' Then we make an outer ring by buffering the inner hull by about one
#' times the expected range. That distance is big enough that, by the 
#' time the field reaches the outer edge, most correlation has died 
#' away—so boundary effects are small. We simplify the outer boundary
#' a little; we don’t need eevery little curve for meshing purposes.

outer_sf   <- st_buffer(inner.hull_sf, 
                        dist = MaxEdge * 5)  
outer_sf_simple <- st_simplify(outer_sf,
                               dTolerance = 0.05 * MaxEdge * 5,
                               preserveTopology = TRUE)

#' Plot inner and outer boundaries
ggplot() +
  geom_sf(data = outer_sf_simple, 
          col = "red",  fill = NA, alpha = 0.5) +
  geom_sf(data = inner.hull_sf,    col = "blue", fill = NA) +
  geom_sf(data = kruger.km, fill = NA, alpha = 0.5) +
  geom_sf(data = schisto.complete.sf, color = "black", size = 0.75) +
  theme_minimal() +
  labs(x = "Easting (km)", y = "Northing (km)")

#' Before we can use the inner and outer boundaries in fm_mesh_2d(), we 
#' need to convert the two sf polygons into their outlines (the lines 
#' that trace around the polygons). The mesher does not take filled 
#' shapes; it only works with these outlines, which it represents 
#' as boundary segments.
inner.hull_segm <- fm_as_segm(inner.hull_sf)
outer.hull_segm <- fm_as_segm(outer_sf_simple)

#' Make the mesh with the new boundary segments
mesh1 <- fm_mesh_2d_inla(
  loc      = Loc,
  boundary = list(inner.hull_segm, outer.hull_segm),
  max.edge = c(1, 5) * MaxEdge,  
  cutoff   = MaxEdge / 5)
fm_crs(mesh1) <- st_crs(schisto.complete.sf)  

#' Size of mesh
mesh1$n 
#' 2852 nodes for 557 observations - decent. We can to set the resolution to finer-scale if we'd like.

#' Ensure the mesh has the same CRS as KNP 
fm_crs(mesh1) <- st_crs(kruger.km)  

#' Plot the study area on the mesh
knp.mesh.plot <- ggplot() +
  theme_minimal() +
  labs(title = "Study area in UTM (km)") +
  geom_fm(data = mesh1) +
  geom_sf(data = schisto.complete.sf,
          alpha = 0.5,
          size = 0.5)  +
  geom_sf(data = kruger.km,
          fill = "transparent", 
          col = "red") 
knp.mesh.plot



# Section 4: Building the Model -----

#* Subsection 4.1: Setting priors -----

#' We need to define the priors for the range of the spatial field, the sigma or
#' u_r(i) of the spatial field. In Bernoulli distributed models, you don't set a
#' prior for the observational noise (the variance associated with the response
#' variable). We will set a PC prior for the RW smoother, though, as we did in v1.

#' Both the north and the south SRFs will use the same priors
#' The range prior is based off of what we found in previous model runs:
#' It is most likely that the range is < 25 km and that the sigma is < 2.

#' After running this model once, it seemed like the range in the south was ~70km
#' and in the north, was about 187km, so we'll go with a prior of 75km, and state
#' that it's most likely that the range is larger than this (P range < 75 is 0.05).
#' The first model run also showed that the st.dev for sigma for the SRFs in 
#' the south in the north were both < 1, so we'll update our prior to reflect
#' that.
Matern <- inla.spde2.pcmatern(mesh1,
                              prior.range = c(75, 0.05),
                              prior.sigma = c(1, 0.05))

#' PC prior for the smoother, initally set to 2, 0.05, but after 1 runthrough,
#' it seemed as though the sd was never higher than 1, so we'll update that
#' for a second runthrough.
PC.rw <- list(prec = list(prior = "pc.prec",
                          param = c(1, 0.05)))
#' This says we think there is only a 5% chance that the sigma_RW is > 2.



#* Subsection 4.2: Defining the smoother -----

#' The `fm_mesh_1d()` function is from the `fmesher` package, and the `bru_mapper()`
#' function is from the `inlabru` package.
#' You can try different values for the number of knots.
mesh.fire <- fm_mesh_1d(seq(min(schisto.complete.sf$weighted_dist_to_fire),
                            max(schisto.complete.sf$weighted_dist_to_fire),
                            length.out = 20)) #' where you set your number of knots
mapper.fire <- bru_mapper(mesh.fire, indexed = FALSE)
#' setting indexed = FALSE ensures that the mapper uses the actual values, rather
#' than indexed (smoother) values, to compute the RW interpolation between knots.
#' Note that the standardized covariates are NOT used in the RW function.


#* Subsection 4.3: Defining the model components -----

#' Step 1: Create a list of all possible covariates you want to model with.
#' These are copy-paste from v1, but with an SRF that varies by region, defined
#' as `SRF.s` and `SRF.n`. We give each region its own intercept, while keeping
#' names of shared components consistent between the two likelihoods.
#' If we wanted, we could also allow any of the other components to differ
#' between regions, if we have reason to believe they may differ biologically.

#' We'll be naming everything with a "4" here, since 3 similar models were ran 
#' in v1.
cmp4 <- ~ Intercept(1) +
  Age(log.age,
      model = "linear") +
  BCS(BCS.Average,
      model = "linear") +
  TWI(topo_wetness_index.std,
      model = "linear") +
  MSI(jan_apr_23_msi.std,
      model = "linear") +
  EVI(jan_apr_23_evi.std,
      model = "linear") +
  soil(clay_percent_30cm.std,
       model = "linear") +
  fire.rw(weighted_dist_to_fire,
          model = "rw2",
          scale.model = TRUE,
          hyper = PC.rw,
          mapper = mapper.fire,
          constr = TRUE) +
  SRF.s(geometry,
      model = Matern) +
  SRF.n(geometry,
        model = Matern)

#' Step 2: Define the likelihoods and options. We'll split the dataset
#' by region `schisto.complete.sf.s` and `schisto.complete.sf.n` so that later on,
#' we can easily extract the residuals based on the "blocks" of south and north 
#' data.

schisto.complete.sf$ID <- seq_len(nrow(schisto.complete.sf))
schisto.complete.sf.s <- subset(schisto.complete.sf, region == "South")
schisto.complete.sf.n <- subset(schisto.complete.sf, region == "North")

#' Create a binomial likelihood for each subset, linking the south-only data to 
#' the south-specific terms (i.e., `SRF.s`), and the same for the north-only 
#' data.

#' first for the south
lik4.s <- bru_obs(result01 ~ Intercept + Age + BCS + TWI + MSI + EVI + soil + 
                      fire.rw + SRF.s,
                   family = "binomial",
                   data = schisto.complete.sf.s)

#' then for the north
lik4.n <- bru_obs(result01 ~ Intercept + Age + BCS + TWI + MSI + EVI + soil +
                      fire.rw + SRF.n,
                   family = "binomial",
                   data = schisto.complete.sf.n)

#' This setup closely mirrors the ideas introduced during the Highland Stats
#' course in terms of the model.matrix with an interaction term. By splitting 
#' the data into south and north and linking each subset to a different set of model 
#' components, we effectively construct separate design-matrix–like mappings for 
#' the two regions. Components with different names (if we had TWI.n and TWI.s) 
#' correspond to different columns in these mappings and are therefore estimated 
#' separately, while components with the same name are shared across likelihoods.

#' Lastly, we'll specify our output options. Because we have two different defined
#' likelihoods, we would include `control.family` as a list with two entries - one
#' per likelihood - for the variation in the response variable (the observational
#' noise, PC.sigma)
Options4 <- list(quantiles =  c(0.025, 0.975),
                control.compute = list(dic = TRUE, waic = TRUE),
                control.inla = list(strategy = "gaussian",
                                    int.strategy = "eb",
                                    diagonal = 1e-6)) 
                #control.family = list(list(hyper = PC.sigma),
                                # list(hyper = PC.sigma))



#* Subsection 4.4: Fitting the model -----

#' We will fit the model with just one call to `bru()`, which will jointly
#' estimate all shared and region-specific components.
mod4 <- bru(components = cmp4, lik4.s, lik4.n, options = Options4)
mod4$summary.fixed[,c(1,3,4)] #' specifically selecting the mean and quantiles




# Section 5: Model Results -----

#* Subsection 5.1: Smoothers -----

#' We'll first plot shared the smoother for weighted distance to fire, then the two 
#' SRFs for each region.

#' #' First, create a grid of distance to fire values between the observed min 
#' and max values.
fire.data.4 <- data.frame(weighted_dist_to_fire = seq(from = min(schisto.complete.sf$weighted_dist_to_fire),
                                                    to = max(schisto.complete.sf$weighted_dist_to_fire),
                                                    length = 100))

#' Now, give this data to the `predict()` function to get outputs for the smoother
pred.fire4 <- predict(mod4, #' tell it which model to use
                      newdata = fire.data.4, #' tell it what data to use in its prediction
                      formula = ~fire.rw, #' tell it what smoother you used
                      n.samples = 1000) #' tell it how many predictions to generate

#' And plot
plot_fire4 <- ggplot(data = pred.fire4,
                     aes(x = weighted_dist_to_fire,
                         y = mean)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_ribbon(aes(ymax = q0.975,
                  ymin = q0.025),
              alpha = 0.2) +
  xlab("Weighted Distance to Fire") + ylab("f(Weighted Distance to Fire)") +
  ggtitle("Model 4") +
  theme(text = element_text(size = 12))
plot_fire4
#' Compared to the other models, this smoother no longer shows an important 
#' effect at mid-distance/intermediate times from fire. This smoother still 
#' indicates an important effect at locations that are further away from fires 
#' that happened >60d ago.

#' To see all 3 model smoother outputs, make sure the models in v1 are loaded
#' before running this line of code:
plot_grid(plot_fire1, plot_fire2, plot_fire4, ncol = 1)


#* Subsection 5.2: SRFs -----

#' Now lets plot the two SRFs to see if they differ from the combined plot.

#' Step 1: Define two new inner and outer hulls to mask the predictions based
#' off of capture locations in each region.

#' Southern region first
loc.s <- cbind(schisto.complete.sf.s$Xkm, schisto.complete.sf.s$Ykm)
inner.hull.s <- fm_nonconvex_hull(loc.s, convex = -0.15, crs = st_crs(schisto.complete.sf.s))
#' Results
plot(inner.hull.s)
points(loc.s)

#' converting to an `sfc` object if we want to mask the SRF to this boundary
InnerHull_sfc.s <- fm_as_sfc(inner.hull.s) 
st_crs(InnerHull_sfc.s) <- st_crs(schisto.complete.sf.s) 

#' Now make the outer hull to use as the SRF boundary
outer_sf.s   <- st_buffer(inner.hull.s, 
                        dist = MaxEdge * 5)  
outer_sf_simple.s <- st_simplify(outer_sf.s,
                               dTolerance = 0.05 * MaxEdge * 5,
                               preserveTopology = TRUE)
OuterHull_sfc.s <- st_geometry(outer_sf_simple.s)

#' Repeat for the north
loc.n <- cbind(schisto.complete.sf.n$Xkm, schisto.complete.sf.n$Ykm)
inner.hull.n <- fm_nonconvex_hull(loc.n, convex = -0.15, crs = st_crs(schisto.complete.sf.n))
#' Results
plot(inner.hull.n)
points(loc.n)

#' converting to an `sfc` object if we want to mask the SRF to this boundary
InnerHull_sfc.n <- fm_as_sfc(inner.hull.n) 
st_crs(InnerHull_sfc.n) <- st_crs(schisto.complete.sf.n) 

#' Now make the outer hull to use as the SRF boundary
outer_sf.n   <- st_buffer(inner.hull.n, 
                          dist = MaxEdge * 5)  
outer_sf_simple.n <- st_simplify(outer_sf.n,
                                 dTolerance = 0.05 * MaxEdge * 5,
                                 preserveTopology = TRUE)
OuterHull_sfc.n <- st_geometry(outer_sf_simple.n)

#' Step 2: Generate the predictions for both regions

NewData.s <- fm_pixels(mesh1, 
                       dims = c(250, 250),
                       mask = OuterHull_sfc.s)

NewData.n <- fm_pixels(mesh1, 
                       dims = c(250, 250),
                       mask = OuterHull_sfc.n)

Pred.s <- predict(object = mod4,
                  newdata = NewData.s,
                  formula = ~SRF.s)

Pred.n <- predict(object = mod4,
                  newdata = NewData.n,
                  formula = ~SRF.n)

#' Define the color palette
ColSc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE))}

#' Step 3: Plot the SRFs

#' First, the south
SRF.plot.s <- ggplot() +
              gg(Pred.s, geom = "tile") +
              geom_sf(data = kruger.km,
                      fill = "transparent", 
                      col = "darkred") +
              gg(schisto.complete.sf.s, size = 0.1) +
              ggtitle("Posterior mean SRF - South KNP") +
              labs(x = "Easting (km)", y = "Northing (km)") + 
              ColSc(Pred.s$mean)
SRF.plot.s

#' Now the north
SRF.plot.n <- ggplot() +
              gg(Pred.n, geom = "tile") +
              geom_sf(data = kruger.km,
                      fill = "transparent", 
                      col = "darkred") +
              gg(schisto.complete.sf.n, size = 0.1) +
              ggtitle("Posterior mean SRF - North KNP") +
              labs(x = "Easting (km)", y = "Northing (km)") + 
              ColSc(Pred.n$mean)
SRF.plot.n

#' Combine the plots
region.combo.SRF <- plot_grid(SRF.plot.s, SRF.plot.n, ncol = 2)
region.combo.SRF

#' Interesting, it seems as though perhaps the hotspot around Skukuza did mask
#' patterns further north in the park - look at the scale of the means in both
#' plots.

#' What's even more interesting is the output when running with updated priors -
#' a HUGE change in what the SRF in the north looks like, not that huge of a 
#' difference in the south.

#* Subsection 5.3: Hyperparameters -----

#' Lets look at the posterior distribution of the range and spatial dependency
#' for both regions.
mod4$summary.hyperpar
#' interestingly, the range is 70km for the south and 187km for the north, and 
#' the 95% credible interval is HUGE in the north.

#' Updating our priors for the second runthrough, the range in the south became
#' ~150km and in the north became ~436km?

#' Step 1: Get the posterior distribution of the range
Range.s <- spde.posterior(mod4, "SRF.s", what = "range")
rangeplot.s <- plot(Range.s) + coord_cartesian(xlim = c(0, 500), ylim = c(0, 0.014)) +
               labs(title = "South KNP Range")

Range.n <- spde.posterior(mod4, "SRF.n", what = "range")
rangeplot.n <- plot(Range.n) + coord_cartesian(xlim = c(0, 500), ylim = c(0, 0.014)) +
               labs(title = "North KNP Range")

region.combo.range <- plot_grid(rangeplot.s, rangeplot.n, ncol = 1)
region.combo.range
#' Interesting - it seems like the peak range is lower for the north, but that 
#' it decreases a bit more gradually than in the south.
#' Updating the priors reaaaallly flattened out the range in both areas, especially
#' in the north!

#' Step 2: Get the imposed spatial dependency by region
MaternCor.s <- spde.posterior(mod4, "SRF.s", what = "matern.correlation")
MaternCor.s$Region <- "South"
pSpatDep.s <- plot(MaternCor.s) + labs(x = "Distance (km)", y = "Correlation") +
              coord_cartesian(ylim = c(0, 1), xlim = c(0, 500)) + 
              labs(title = "South KNP")
pSpatDep.s

MaternCor.n <- spde.posterior(mod4, "SRF.n", what = "matern.correlation")
MaternCor.n$Region <- "North"
pSpatDep.n <- plot(MaternCor.n) + labs(x = "Distance (km)", y = "Correlation") +
              coord_cartesian(ylim = c(0, 1), xlim = c(0, 500)) + 
              labs(title = "North KNP")
pSpatDep.n

#' Combine the plots
MaternCor.combo <- rbind(MaternCor.s, MaternCor.n)
MaternCor.combo <- as.data.frame(MaternCor.combo)
MaternCor.combo$Region <- as.factor(MaternCor.combo$Region)

plot.combo.SpatDep <- ggplot(data = MaternCor.combo,
                             aes(x = x,
                                 y = q0.5,
                                 group = Region,
                                 col = Region)) +
                      geom_line() + 
                      labs(x = "Distance (km)", y = "Correlation") +
                      coord_cartesian(ylim = c(0,1), xlim = c(0,500))
plot.combo.SpatDep
#' Again, you can see that the spatial dependency in the north decays more slowly
#' with distance, indicating a stronger and more persistent spatial structure
#' than in the south. In the south, spatial correlation drops rapidly - it is 
#' more local - and is closer to 0 at shorter distances.
#' We could include the code that defines where the correlation drops below
#' certain thresholds to get actual numbers and see differences more clearly
#' between regions, if we'd like.

#' With the updated priors for the second runthrough, the decay became more gradual
#' for both regions, but especially for the north!

#* Subsection 5.4: Validation and fit -----
#' Make sure the models in v1 are loaded before running this next bit.

#' Step 1: Extract DIC and WAIC and add to the existing values for the v1 models
DICm4 <- mod4$dic$dic
WAICm4 <- mod4$waic$waic

DICs <- c(DICs, DICm4)
WAICs <- c(WAICs, WAICm4)

prevOutput <- data.frame(Models = c("GAM", "GAM + SRF", "GLM + Age", "GAM + SRF by Region v1", 
                                    "GAM + SRF by Region v2"),
                  DIC = DICs,
                  WAIC = WAICs)
colnames(prevOutput) <- c("      Model",
                          "        DIC",
                          "       WAIC")
prevOutput
#' It seems as though including a different SRF by region doesn't really improve
#' the model in terms of DIC/WAIC. The lowest WAIC, interestingly, is the model
#' with the GAM and single SRF.
#' SRF by region *is* however, an improvement on the model with no SRF at all.
#' The second runthrough with updated priors based off of the first runthrough 
#' has a higher WAIC, but is still lower than the GAM without SRF.

#' Step 2: Model validation
#' Even though this model isn't really an improvement over models 2 and 3, it's
#' still important to check model fit.

#' By splitting up the data in the north and south, and applying a multiple 
#' likelihood model, the fitted values mod4$summary.fitted.values are now in 
#' blocks. The first block is for the south (we typed them first when 
#' defining the likelihood) and the second block is for the north. We can 
#' extract these two blocks.
ns <- nrow(schisto.complete.sf.s)
nn <- nrow(schisto.complete.sf.n)
Fit.s <- mod4$summary.fitted.values[1:ns, "mean"]
Fit.n <- mod4$summary.fitted.values[(ns + 1):(ns + nn), "mean"]

#' We then add these fitted values at the appropriate rows using the ID 
#' column we added earlier in this exercise.
schisto.complete.sf$Fit4 <- NA
schisto.complete.sf$Fit4[match(schisto.complete.sf.s$ID, schisto.complete.sf$ID)] <- Fit.s
schisto.complete.sf$Fit4[match(schisto.complete.sf.n$ID, schisto.complete.sf$ID)] <- Fit.n
#' Now we have the fitted values in the right order inside the `schisto.complete.sf` 
#' object. And that means that we can get our scaled quantile residuals.

#' Start a loop simulating 1000 Bernoulli distributed observations for each 
#' given fitted value, Fit4, per likelihood. We'll use the `rbinom()` function 
#' from the `stats` package.
N4 <- nrow(schisto.complete.sf)
Ysim4 <- matrix(nrow = N4, ncol = 1000)
for (i in 1:1000){
  Ysim4[,i] <- rbinom(n = N4,
                      size = 1,
                      prob = schisto.complete.sf$Fit4)
}

#' Now pass this simulated data through the `DHARMa` package, which will calculate
#' scaled quantile residuals for you.
E4sqr <- createDHARMa(simulatedResponse = Ysim4,
                      observedResponse = schisto.complete.sf$result01,
                      fittedPredictedResponse = schisto.complete.sf$Fit4,
                      integerResponse = TRUE)

#' If you were running a Poisson model, here's where you would check for
#' overdispersion.

#' Step 2: Check for residual uniformity
#' The `DHARMa` package will also do this for you!
par(mfrom = c(1,1), mar = c(5,5,2,2))
plotQQunif(E4sqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)
#' Seems ok - we may be overfitting?

#' Step 3: Plot the residuals vs the covariates using `DHARMa`
plotResiduals(E4sqr, form = schisto.complete.sf$log.age) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$BCS.Average) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$topo_wetness_index.std) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$clay_percent_30cm.std) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$jan_apr_23_evi.std) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$jan_apr_23_msi.std) #' Ok
plotResiduals(E4sqr, form = schisto.complete.sf$weighted_dist_to_fire) #' Ok
#' all seem ok here.

#' Step 4 (Optional): Plot residuals vs space
schisto.complete.sf$E4 <- residuals(E4sqr)

#' Define the color and point size for the plot
#' Remember that scaled quantile residuals are centered on 0.5, while Pearson 
#' residuals are centered around 0.
schisto.complete.sf$MyCol <- ifelse(schisto.complete.sf$E4 >= 0.5, "red", "blue")
schisto.complete.sf$MySize <- rescale(abs(schisto.complete.sf$E4), to = c(0,3))

#' Plot it
mod4residspat <- ggplot(data = kruger.km) +
  geom_sf(fill = "transparent") +
  geom_sf(data = schisto.complete.sf,
          aes(col = MyCol,
              size = MySize),
          alpha = 0.3) +
  scale_color_identity() +
  scale_size_continuous(range = c(1,3)) +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Residuals Mod 4")
mod4residspat
#' This looks the same as the others
plot_grid(mod1residspat, mod2residspat, mod3residspat, mod4residspat, ncol = 2)
#' Very little, if any, difference between mod 3 output

#' Step 5: Create a semi-variogram
#' First, make a dataframe with the residuals and capture locations
mod4data <- data.frame(E4 = schisto.complete.sf$E4,
                       Xkm = schisto.complete.sf$Xkm,
                       Ykm = schisto.complete.sf$Ykm)

#' Convert to an `sf` object
mod4data_sf <- st_as_sf(x = mod4data,
                        coords = c("Xkm", "Ykm"),
                        crs = NA) #' CRS doesn't matter for this, so no need to specify

#' Call the `variogram()` function from the package `gstat`
V4 <- variogram(E4 ~1,
                data = mod4data_sf,
                cutoff = 100,
                cressie = TRUE)

#' And plot
variogrammod4 <- ggplot(data = V4, aes(x = dist, y = gamma)) +
  geom_point() +
  geom_smooth(se = FALSE, span = 1) +
  labs(x = "Distance (km)", y = "Semi-variogram Mod 4") +
  coord_cartesian(ylim = c(0.083, 0.093)) +
  theme(text = element_text(size = 15))
variogrammod4

plot_grid(variogrammod1, variogrammod2, variogrammod3, variogrammod4, ncol = 2)
#' This variogram is a bit less steep than that for mod 3.

#* Subsection 5.5: Fixed effects -----

#' First runthrough:
#' Let's take a look at the posterior mean and 95% credible intervals for the
#' intercept and covariates
mod4$summary.fixed

#' This means that, for our model in the south we have:
#' logit(Pi_r) = -4.822 - 0.049 x Age + 0.568 x BCS - 0.128 x TWI + 0.116 x MSI +
#'                  0.028 x EVI - 0.032 x Soil + f(fire) + SRF.s

#' And in the north:
#' logit(Pi_r) = -4.822 - 0.049 x Age + 0.568 x BCS - 0.128 x TWI + 0.116 x MSI +
#'                  0.028 x EVI - 0.032 x Soil + f(fire) + SRF.n

#' Second runthrough:
#' In the south:
#' logit(Pi_r) = -4.262 - 0.035 x Age + 0.486 x BCS - 0.131 x TWI + 0.080 x MSI +
#'                  0.022 x EVI - 0.020 x Soil + f(fire) + SRF.s
#' Same in the north, just with `SRF.n` 

#' Updating our priors changed the intercept and the coefficients a bit, mostly
#' for MSI, EVI, and soil - perhaps the SRF "ate" up more of the explanation in 
#' variance.