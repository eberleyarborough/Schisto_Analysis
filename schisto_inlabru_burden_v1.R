#########
#' This R script is a variation of the code found in `schisto_inlabru_v1.Rmd`,
#' and seeks to run the same general analyses, but with Schistosoma spp. 
#' *burden,* rather than postiive/negative disease status. 
#########


# Section 1: Setup -----
#' We need the following packages:
library(lattice)
library(ggplot2)
library(sf)
library(gstat)
library(Matrix)
library(INLA)
library(fmesher)
library(inlabru)
library(DHARMa)
library(dplyr)
library(glmmTMB)
library(rnaturalearth)
library(scales)
library(cowplot)
library(tidyr)
source("/Users/eberleyarborough/Documents/Coding Help/INLA tutorials/highland_stats_INLA_course/AllRCodeInlabru/HighstatLibV14.R") #' Highland Statistics support file

#' Set the working directory and import the data (should be the project directory).
setwd("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/Schisto/Schisto_Analysis")


# Section 2: Importing and Transforming Data -----

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

#' We now transform the spatial object to the kilometre-based coordinate 
#' system and extract the coordinates as numeric variables. These will 
#' be used later when building spatial models.
schisto2.sf <- st_transform(schisto.sf, CRS.km)
Coords <- st_coordinates(schisto2.sf)
schisto2.sf$Xkm <- Coords[, 1]
schisto2.sf$Ykm <- Coords[, 2]

#' Dropping unnecessary columns from the schisto data frame
schisto2.sf <- schisto2.sf %>%
  select("animal.ID", "age.years", "age.class", "Sex", "BCS.Average", "CAA.Estimate", "Result", "Xkm", "Ykm", "geometry")

#' Step 2: Import the environmental covariates
#' Read in the .csv file with the environmental covariates. NOTE: This file
#' contains data for all capture locations, not just those that are schisto positive.
env_data <- read.csv("/Users/eberleyarborough/Documents/2023 Buffalo XS Disease Data/CSVs/cap_and_env_data.csv")
env_data <- env_data %>%
  select("animal.ID", "capture.date", "capture.day", "capture.group.by.day", "approx.herd.size", "weighted_dist_to_fire", "topo_wetness_index", "clay_percent_30cm", "jan_apr_23_evi", "jan_apr_23_msi")

#' Join to schisto data
schisto.complete.sf <- left_join(schisto2.sf, env_data, by = "animal.ID")

#' Step 3: Transform Data
#' The cutoff for a positive schisto test on CAA is a titer over 30, so we will
#' drop rows where the CAA titer is =< 20.
schisto.pos.sf <- schisto.complete.sf %>%
  dplyr::filter(CAA.Estimate >= 30)
#' We're only left with 101 observations.

#' And convert Sex, age class, and capture group by day to factors
schisto.pos.sf$fage.class <- factor(schisto.pos.sf$age.class)
schisto.pos.sf$fSex <- factor(schisto.pos.sf$Sex)
schisto.pos.sf$fcapgroup <- factor(schisto.pos.sf$capture.group.by.day)

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
  geom_sf(data = schisto.pos.sf) +
  theme_minimal()

#' Confirm that the data and the map share the same CRS: 
st_crs(kruger.km) == st_crs(schisto.pos.sf)
#' Output should be "TRUE"


# Section 3: Building the Mesh -----

#' Here's where we define our mesh components and construct the mesh that we'll
#' use later on for our spatial random field (SRF).

#' Step 1: Making the mesh -----
#' To get a good idea of mesh settings for the SRF, we first need to get a sense what the distances are between the 
#' sampling locations. 
Loc <- cbind(schisto.pos.sf$Xkm, schisto.pos.sf$Ykm)
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
#' We'll use the same resolution as we used for our second runthrough of the
#' Bernoulli schisto model.
MeshResolution <- 15               #' in km.
MaxEdge  <- MeshResolution / 5    #' good rule of thumb is to have max edge set to 1/5 the resolution

#' #' fm_nonconvex_hull makes a non-convex area around the sampling
#' locations - the inner boundary. You can control how close the blue line is 
#' to the sites with the 'convex' argument. 
inner.hull <- fm_nonconvex_hull(Loc, convex = -0.05, crs = st_crs(schisto.pos.sf))
#' Results
plot(inner.hull)
points(Loc)

#' Convert the inner, nonconvex hull to an sf object
inner.hull_sf <- st_as_sf(inner.hull)

#' Then we make an outer ring by buffering the inner hull by about one
#' times the expected range. That distance is big enough that, by the 
#' time the field reaches the outer edge, most correlation has died 
#' away—so boundary effects are small. We simplify the outer boundary
#' a little; we don’t need every little curve for meshing purposes.

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
  geom_sf(data = schisto.pos.sf, color = "black", size = 0.75) +
  theme_minimal() +
  labs(x = "Easting (km)", y = "Northing (km)")
#' These looks nearly identical to the boundaries for the Bernoulli model.

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
fm_crs(mesh1) <- st_crs(schisto.pos.sf)  

#' Size of mesh
mesh1$n 
#' 6173 nodes for 101 observations - decent, bordering on too fine.

#' Ensure the mesh has the same CRS as KNP 
fm_crs(mesh1) <- st_crs(kruger.km)  

#' Plot the study area on the mesh
knp.mesh.plot <- ggplot() +
  theme_minimal() +
  labs(title = "Study area in UTM (km)") +
  geom_fm(data = mesh1) +
  geom_sf(data = schisto.pos.sf,
          alpha = 0.5,
          size = 0.5)  +
  geom_sf(data = kruger.km,
          fill = "transparent", 
          col = "red") 
knp.mesh.plot


# Section 4: Data Exploration -----
#' Now we'll explore the data a bit more to look at any underlying relationships, 
#' nonlinear patterns, etc.

#' Step 1: Determine if there are any missing values
MyVar <- c("animal.ID", "age.years", "age.class", "Sex", "BCS.Average", "CAA.Estimate",
           "capture.date", "capture.day", "capture.group.by.day", "approx.herd.size",
           "weighted_dist_to_fire", "topo_wetness_index", "clay_percent_30cm",
           "jan_apr_23_evi", "jan_apr_23_msi")
colSums(is.na(schisto.pos.sf[,MyVar]))
#' Looks like missing values in herd size, EVI, and MSI (4, 6, and 6, respectively).

#' We'll drop NA rows
#' Now we'll remove any NA columns
schisto.pos.sf <- schisto.pos.sf %>%
  drop_na()
#' Left with 91 observations

#' For some analyses, it is convenient to work with a standard data frame 
#' rather than an sf object. We therefore drop the geometry column and 
#' store the result separately.
schisto.pos.df <- st_drop_geometry(schisto.pos.sf)

#' It's also a good idea to standardize continuous covariates - this can only be done after NA is dropped.
#' The `MyStd()` function is a custom function in the Highland Statistics 
#' script.
schisto.pos.sf$weighted_dist_to_fire.std <- MyStd(schisto.pos.sf$weighted_dist_to_fire)
schisto.pos.sf$topo_wetness_index.std <- MyStd(schisto.pos.sf$topo_wetness_index)
schisto.pos.sf$clay_percent_30cm.std <- MyStd(schisto.pos.sf$clay_percent_30cm)
schisto.pos.sf$jan_apr_23_evi.std <- MyStd(schisto.pos.sf$jan_apr_23_evi)
schisto.pos.sf$jan_apr_23_msi.std <- MyStd(schisto.pos.sf$jan_apr_23_msi)

#' Step 2: Look the depth of the data - do we have enough observations with the
#' given covariates?

#' How many observations per date?
table(schisto.pos.sf$capture.date)
#' Anywhere from 1 to 14

#' How many groups per day?
table(schisto.pos.sf$capture.group.by.day)
#' as much as 4 groups per day - consider using as a random effect?

#' What's the range of herd sizes?
hist(schisto.pos.sf$approx.herd.size)
#' looks like 50 - 600 animals per herd - consider using as a random effect?

#' Are herd size and capture group per day related?
plot(x = schisto.pos.sf$capture.group.by.day,
     y = schisto.pos.sf$approx.herd.size,
     xlab = "Capture Group by Day",
     ylab = "Approx. Herd Size")
#' fairly even spread of herd sizes across capture groups - pick one to use
#' as a random effect, not both

#' How many per age class?
table(schisto.pos.sf$age.class)
#' 15 calves, 58 subadults, 16 adults, 2 old adults - would need to combine
#' adults and old adults if wanting to use this as a covariate or as an 
#' interaction term

#' What is the age range?
hist(schisto.pos.sf$age.years)
range(schisto.pos.sf$age.years)
#' as young as 0.66 to 20 years

#' What are the differences in sex?
table(schisto.pos.sf$Sex) #' enough observations between sexes, though just barely

#' Step 3: Check for outliers - the `MyDotplot.ggp2()` function comes from the
#' Highland Statistics custom script.
MyDotplot.ggp2(schisto.pos.df, MyVar, Ncol = 3, TextSize = 10,
               PointSize = 0.5, DropLabels = FALSE)
#' No obvious outliers, but possible collinearity or non-linear patterns.

#' Step 4: Assess for collinearity
#' Make a pairplot - the `Mypairs()` function comes from the Highland Statistics
#' custom script.
MyVar <- c("Xkm", "Ykm", "age.years", "BCS.Average", "clay_percent_30cm",
           "jan_apr_23_evi", "jan_apr_23_msi", "topo_wetness_index",
           "weighted_dist_to_fire")
Mypairs(schisto.pos.df[,MyVar])
#' There might be some collinearity between Xkm and/or Ykm and percentage of clay
#' in the soil, as well as EVI and MSI, which makes sense, as these are spatial
#' covariates.

#' We only should be nervous if a Pearson correlation is larger than 0.8; however,
#' that doesn't mean we should ignore moderately high (~0.6) values, like
#' those for X/Ykm and clay percentage.

#' If we had temporal (longitudinal) data, here's where we would check to see how/if
#' covariates change over time by plotting every covariate vs. year.

#' Do covariates differ by clay percentage in the soil? We anticipate there 
#' may be a relationship, as the clay content influences water retention
#' on the landscape (MSI) and vegetation growth (EVI).
#' First, we need to put the data into "long format" for easier plotting.
GatherTheseCols <- c("jan_apr_23_msi", "jan_apr_23_evi", "topo_wetness_index",
                     "weighted_dist_to_fire", "age.years")
schist_long <- gather(data = schisto.pos.sf,
                      key = "ID",
                      value = "AllX",
                      all_of(GatherTheseCols),
                      factor_key = TRUE)
dim(schist_long)

#' Check to see how much variation in EVI/MSI is explained by clay percentage
#' with SLR. There wasn't much of a relationship in the Bernoulli model, but
#' just a sanity check here.
test_msi <- lm(jan_apr_23_msi ~ clay_percent_30cm, data = schisto.pos.sf)
summary(test_msi)
#' Less than 1%

test_evi <- lm(jan_apr_23_evi ~ clay_percent_30cm, data = schisto.pos.sf)
summary(test_evi)
#' Less than 1%

#' Step 5: Continuous covariates vs. disease data
#' Plot each continuous covariate against the CAA titer data. Note: the GAM
#' smoother function family should match the type of model you're planning
#' on running (e.g., negative binomial, gaussian, etc.)
ggplot(data = schist_long, aes(x=AllX, y = CAA.Estimate)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Covariate", y = "Disease Status") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 5)) +
  geom_smooth(method = "gam",
              method.args = list(family = "gaussian"),
              se = FALSE) +
  facet_wrap(.~ID, scale = 'free_x', ncol = 2)
#' Not sure how relevant this is if we're not running a non-gaussian model.

#' Step 6: Assess for dependency
#' Here are the spatial patterns of the capture locations (using base R)
xyplot(Ykm ~ Xkm,
       aspect = "iso",
       col = 1,
       pch = 16,
       data = schisto.pos.sf)

xyplot(Ykm ~ Xkm | fcapgroup,
       aspect = "iso",
       col = 1,
       pch = 16,
       data = schisto.pos.sf)
#' It looks like extra captures (like 1a 1b, 2a 2b) happened in the north,
#' and up to 4 captures per day happened in the south.
#' Most capture days had only 2 captures per day. Just something to think about
#' if we wanted to include capture group by day in the analyses in any way. We 
#' would probably remove 1a/b and 2a/b since there's no real variation there, and
#' we would just call them "1" and "2."

#' Here's the same kind of plot as a ggplot
ggplot(data = kruger.km) + #' boundary
  geom_sf(fill = "transparent") + #' no fill
  geom_sf(data = schisto.pos.sf, #' capture locations
          size = 0.75, #' size of the points
          alpha = 0.5) + #' transparency of the points
  theme_minimal() +
  labs(title = "Study Area in UTM (km)")

#' And split up by capture group
ggplot(data = kruger.km) +
  geom_sf(fill = "transparent") +
  geom_sf(data = schisto.pos.sf,
          size = 0.75,
          alpha = 0.5) +
  theme_minimal() +
  labs(title = "Study Area in UTM (km)") +
  facet_wrap(~fcapgroup)

#' To conclude about our data exploration:
#' - We have spatial patterns in some of our covariates, which isn't surprising.
#' - Herd size or capture group could be used as random effects; however, if using
#' capture group by day, we should collapse 1a/b and 2a/b into 1 and 2.
#' - If we want to use age class, we should collapse "adult" and "old adult" since
#' there are only 3 "old adult" animals.
#' - It doesn't look like we have any dependency (other than potential spatial
#' dependency).


# Section 5: Model Formulas -----

#' Let CAA_i be the CAA titer value at site i.
#' We assume a normal distribution for CAA_i:

#'   CAA_i ~ N(mu_i, sigma^2)
#'   E[CAA_i]   = mu_i
#'   var[CAA_i] = sigma^2

#' **Model 1**
#' We model the mean mu_i as a function of the host demographics, environmental
#' covariates, and spatial dependency.

#' mu_i = Intercept + age.years + BCS.Average + TWI_i + EVI_i + MSI_i + soil_i +
#'        firedistance_i + u_i

#' u ~ N(0, Omega)

#' **Model 2**
#' To assess whether we really need the spatial term u_i, we also execute 
#' a model without spatial dependency. That is this model:

#' mu_i = Intercept + age.years + BCS.Average + TWI_i + EVI_i + MSI_i + soil_i +
#'        firedistance_i


# Section 6: Building the Models

#' We'll first fit the model *without* spatial dependency and use that model as 
#' the reference model.

#' In `inlabru()`, there are always 3 steps before fitting a model:
#'    1. Defining the components - what is going into the model (the "suitcase")
#'    2. Defining the likelihood/model formula
#'    3. Defining output options (decide what numerical output you'd like)

#' Keep in mind that `inlabru()` treats factors as random effects in the output.
#' Also note that when defining components, you can name whatever you'd like 
#' *outside* the parenthesis for that covariate, but what is *inside* the 
#' parentheses has to match the name in the data frame.

#' Step 1: Assign priors
#' Setting the PC priors for range and sigma of the SRF. This uses the mesh 
#' we made earlier and the `INLA` package.
#' We'll state that it is most likely that the range is smaller than 50km, given
#' what we know from our Bernoulli model.

#' For a Gaussian model with an identity link, we can use the sd of the response
#' variable as a lower limit of the sigma prior.
sd(schisto.pos.sf$CAA.Estimate)

#' Wang et al (2018) use 1.5 (or 3?) times the sd for the lower 
#' limit in the PC prior.
#' Let's go with twice the estimated SD. This would mean that we think that it
#' is unlikely that the sigma is larger than  7000.
Matern <- inla.spde2.pcmatern(mesh1,
                              prior.range = c(50, 0.95),
                              prior.sigma = c(7000, 0.05))

#' Then there is also the observation noise. That is the sigma in:
#' CAA_i ~ N(mu_i, sigma^2). We will also use a PC prior:
#' P(sigma > 1) = 0.05 
PCPrior.Sigma <- list(prec = list(prior = "pc.prec", 
                                  param = c(1, 0.05)))

#' Step 2: Define components
#' Model 1: LM **without** spatial dependency
#' You can put anything you think you *might* need or want in your components
#' list, but you don't have to include them all in your likelihood.
cmp1 <- ~ Intercept(1) +
  Age(age.years,
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
  fire(weighted_dist_to_fire.std,
          model = "linear")

#' If including a factor, you would write that as:
#' fcap.group(fcapgroup,
#'            model = "factor_contrast")

#' Model 2: LM **with** spatial dependency
#' This is copy-paste from above, just with the SRF added at the end
cmp2 <- ~ Intercept(1) +
  Age(age.years,
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
  fire(weighted_dist_to_fire.std,
          model = "linear") +
  SRF(geometry,
      model = Matern)

#' Step 3: Define likelihood
#' Now we'll define our likelihood, which is the same as specifying the model 
#' formula. Here you can pick any or all of the components we defined in the 
#' previous step. Just keep in mind that here you need to use the covariate 
#' names that we defined above *outside* of the parentheses.

#' Model 1: LM **without** spatial dependency
#' Note that the `bru_obs()` command comes from the `inlabru` package
lik1 <- bru_obs(CAA.Estimate ~ Intercept + Age + BCS + TWI + MSI + EVI + soil + 
                  fire,
                options = list(control.family = list(hyper = PCPrior.Sigma)),
                family = "gaussian",
                data = schisto.pos.sf)

#' Model 2: LM **with** spatial dependency
lik2 <- bru_obs(CAA.Estimate ~ Intercept + Age + BCS + TWI + MSI + EVI + soil + 
                  fire + SRF,
                options = list(control.family = list(hyper = PCPrior.Sigma)),
                family = "gaussian",
                data = schisto.pos.sf)

#' Step 4: Specify output options
#' Now quickly define the output options for the models, which will all be the 
#' same. We would like the 95% credible intervals for the covariates, the DIC, 
#' and the WAIC for the models. Specifying this here means you won't get all the 
#' extra output info you don't want.
#' This will be the same across models
Options <- list(quantiles =  c(0.025, 0.975),
                control.compute = list(dic = TRUE, waic = TRUE))


# Section 7: Fitting the Models -----

#' Model 1: LM **without** spatial dependency
#' Step 1: Fit the model
mod1 <- bru(components = cmp1, lik1, options = Options)

#' Step 2: Look at initial results
mod1$summary.fixed
mod1$summary.hyperpar
#' for random effects, you would run the command:
#' mod1$summary.random$...
#' The precision for the observation noise sigma is 6.78

#' Model 2: BLM **with** spatial dependency
#' Step 1: Fit the model
mod2 <- bru(components = cmp2, lik2, options = Options)

#' Step 2: Look at initial results
mod2$summary.fixed
mod2$summary.hyperpar
#' for random effects, you would run the command:
#' mod1$summary.random$...
#' The precision for the observation noise sigma is still 6.78
#' The range for the SRF is 26.58, and the sd is 64.67

#' Including the SRF changed a lot about what was deemed "important" in terms
#' of covariate credible intervals.


# Section 8: Model Validation -----

#' Here we'll look at scaled quantile residuals, uniformity in residuals, and 
#' check spatial dependency for each model.

#' Model 1: LM **without** spatial dependency
#' Step 1: Get scaled quantile residuals
#' To do this, we first need to get the fitted values from the model
N <- nrow(schisto.pos.sf)
schisto.pos.sf$Pi.mod1 <- mod1$summary.fitted.values[1:N, "mean"]

#' And we need the sigma for the "family noise" or the residual standard
#' deviation. Remember that we can't simply back transform the precision (tau)
#' to sigma directly.
tau1     <- mod1$marginals.hyperpar$`Precision for the Gaussian observations`
MySqrt  <- function(x) { 1 / sqrt(x) }
Sigma1   <- inla.emarginal(MySqrt, tau1)
Sigma1
#' 3872.235, close to what we estimated it would be

#' Start a loop simulating 1000 gaussian distributed observations in the 
#' response variable, for each  given fitted value, Pi.mod1. We'll use the 
#' `rnorm()` function from the `stats` package. Note: gaussian simulations
#' require the mean fitted values and residual standard deviation.
Ysim1 <- matrix(nrow = N, ncol = 1000)
for (i in 1:1000){
  Ysim1[,i] <- rnorm(n = N,
                      mean = schisto.pos.sf$Pi.mod1,
                     sd = Sigma1)
}

#' Now pass this simulated data through the `DHARMa` package, which will calculate
#' scaled quantile residuals for you.
E1sqr <- createDHARMa(simulatedResponse = Ysim1,
                      observedResponse = schisto.pos.sf$CAA.Estimate,
                      fittedPredictedResponse = schisto.pos.sf$Pi.mod1,
                      integerResponse = FALSE) #' is "TRUE" for Bernoulli and poisson data

#' If you were running a Poisson model, here's where you would check for
#' overdispersion.

#' Step 2: Check for residual uniformity
#' The `DHARMa` package will also do this for you!
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1sqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)
#' Woof - not ok!

#' Step 3: Plot the residuals vs the covariates using `DHARMa`
plotResiduals(E1sqr, form = schisto.pos.sf$age.years) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$BCS.Average) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$topo_wetness_index.std) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$clay_percent_30cm.std) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$jan_apr_23_evi.std) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$jan_apr_23_msi.std) #' Nope
plotResiduals(E1sqr, form = schisto.pos.sf$weighted_dist_to_fire) #' Nope
#' Lots of problems here!!

#' Lets also look at the deviance residuals (CAA_i - Fitted_i) / sigma using
#' the `residuals()` function from the `stats` package.
schisto.pos.sf$deviance.mod1 <- residuals(mod1)$deviance.residuals

#' We can also calculate the ordinary residuals, which are just CAA_i - Fitted_i
#' for a Gaussian model.
schisto.pos.sf$resid.mod1 <-  schisto.pos.sf$CAA.Estimate - schisto.pos.sf$Pi.mod1 #' obs - fitted

#' Now we can check for homogeneity of variance plotting the ordinary residuals
#' vs. fitted values
pResidvsFit <- ggplot(data = schisto.pos.sf, aes(y = resid.mod1, x = Pi.mod1)) + 
  geom_point(shape = 1, size = 1) +
  labs(x = "Fitted values", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsFit
#' Problems here

#' And we can check the normality of the ordinary residuals
pResidNorm <- ggplot(data = schisto.pos.sf, aes(x = resid.mod1)) + 
  geom_histogram(alpha = 0.5, color = "darkblue", fill = "lightblue") +
  labs(y = "Histogram", x = "Residuals") +
  theme(text = element_text(size=15))
pResidNorm
#' Not normal in the slightest - need to consider transforming the CAA titers.

#' Based on the scaled quantile residuals, we can probably assume that the plots
#' of residuals vs covariates will have issues, but let's just see one or two.

#' Plot residuals vs. age.years
pResidvsAge <- ggplot(data = schisto.pos.sf, aes(y = resid.mod1, x = age.years)) + 
  geom_point(shape = 1, size = 1) +
  geom_smooth() +
  labs(x = "Age (years)", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsAge
#" Not ideal

#' Plot residuals vs. EVI
pResidvsEVI <- ggplot(data = schisto.pos.sf, aes(y = resid.mod1, x = jan_apr_23_evi.std)) + 
  geom_point(shape = 1, size = 1) +
  geom_smooth() +
  labs(x = "EVI (std)", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsEVI
#' Not too bad, but not great either.

#' Step 4: Check for spatial dependency
#' There are 3 ways we can do this:
#'  1. Plot residuals vs spatial locations and look for patterns
#'  2. Apply a Moran's I test on the residuals
#'  3. Make a variogram of the residuals
#' We'll do options 1 and 3 here.

#' Option 1: Plot residuals vs spatial locations

#' Define the color and point size for the plot
#' Remember that scaled quantile residuals are centered on 0.5, while Pearson 
#' residuals are centered around 0.
schisto.pos.sf$MyCol <- ifelse(schisto.pos.sf$resid.mod1 >= 0, "red", "blue")
schisto.pos.sf$MySize <- rescale(abs(schisto.pos.sf$resid.mod1), to = c(0,10))

#' Plot it
ggplot(data = kruger.km) +
  geom_sf(fill = "transparent") +
  geom_sf(data = schisto.pos.sf,
          aes(col = MyCol,
              size = MySize),
          alpha = 0.3) +
  scale_color_identity() +
  scale_size_continuous(range = c(1,3)) +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Residuals")
#' Hmm, maybe a bit of a spatial pattern?

#' Option 3: Make a variogram
#' First, make a dataframe with the residuals and capture locations
mod1data <- data.frame(E1 = schisto.pos.sf$resid.mod1,
                       Xkm = schisto.pos.sf$Xkm,
                       Ykm = schisto.pos.sf$Ykm)

#' Convert to an `sf` object
mod1data_sf <- st_as_sf(x = mod1data,
                        coords = c("Xkm", "Ykm"),
                        crs = NA) #' CRS doesn't matter for this, so no need to specify

#' Call the `variogram()` function from the package `gstat`
V1 <- variogram(E1 ~1,
                data = mod1data_sf,
                cutoff = 100,
                cressie = TRUE)

#' And plot
ggplot(data = V1, aes(x = dist, y = gamma)) +
  geom_point() +
  geom_smooth(se = FALSE, span = 1) +
  labs(x = "Distance (km)", y = "Semi-variogram") +
  theme(text = element_text(size = 15))
#' Negative spatial correlation peak at ~60km?


#' Model 2: LM **with** spatial dependency
#' Step 1: Get scaled quantile residuals
#' To do this, we first need to get the fitted values from the model
N <- nrow(schisto.pos.sf)
schisto.pos.sf$Pi.mod2 <- mod2$summary.fitted.values[1:N, "mean"]

#' And we need the sigma for the "family noise" or the residual standard
#' deviation. Remember that we can't simply back transform the precision (tau)
#' to sigma directly.
tau2     <- mod2$marginals.hyperpar$`Precision for the Gaussian observations`
MySqrt  <- function(x) { 1 / sqrt(x) }
Sigma2   <- inla.emarginal(MySqrt, tau2)
Sigma2
#' Only 19! So much was taken up by adding the SRF!

#' Start a loop simulating 1000 gaussian distributed observations in the 
#' response variable, for each  given fitted value, Pi.mod1. We'll use the 
#' `rnorm()` function from the `stats` package. Note: gaussian simulations
#' require the mean fitted values and residual standard deviation.
Ysim2 <- matrix(nrow = N, ncol = 1000)
for (i in 1:1000){
  Ysim2[,i] <- rnorm(n = N,
                     mean = schisto.pos.sf$Pi.mod2,
                     sd = Sigma2)
}

#' Now pass this simulated data through the `DHARMa` package, which will calculate
#' scaled quantile residuals for you.
E2sqr <- createDHARMa(simulatedResponse = Ysim2,
                      observedResponse = schisto.pos.sf$CAA.Estimate,
                      fittedPredictedResponse = schisto.pos.sf$Pi.mod2,
                      integerResponse = FALSE) #' is "TRUE" for Bernoulli and poisson data

#' If you were running a Poisson model, here's where you would check for
#' overdispersion.

#' Step 2: Check for residual uniformity
#' The `DHARMa` package will also do this for you!
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E2sqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)
#' this is even WORSE compared to model 1!

#' Step 3: Plot the residuals vs the covariates using `DHARMa`
plotResiduals(E2sqr, form = schisto.pos.sf$age.years) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$BCS.Average) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$topo_wetness_index.std) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$clay_percent_30cm.std) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$jan_apr_23_evi.std) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$jan_apr_23_msi.std) #' Nope
plotResiduals(E2sqr, form = schisto.pos.sf$weighted_dist_to_fire) #' Nope
#' Lots of problems here!!

#' Lets also look at the deviance residuals (CAA_i - Fitted_i) / sigma using
#' the `residuals()` function from the `stats` package.
schisto.pos.sf$deviance.mod2 <- residuals(mod2)$deviance.residuals

#' We can also calculate the ordinary residuals, which are just CAA_i - Fitted_i
#' for a Gaussian model.
schisto.pos.sf$resid.mod2 <-  schisto.pos.sf$CAA.Estimate - schisto.pos.sf$Pi.mod2 #' obs - fitted

#' Now we can check for homogeneity of variance plotting the ordinary residuals
#' vs. fitted values
pResidvsFit <- ggplot(data = schisto.pos.sf, aes(y = resid.mod2, x = Pi.mod2)) + 
  geom_point(shape = 1, size = 1) +
  labs(x = "Fitted values", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsFit
#' Problems here, but does look a bit better than model 1.

#' And we can check the normality of the ordinary residuals
pResidNorm <- ggplot(data = schisto.pos.sf, aes(x = resid.mod2)) + 
  geom_histogram(alpha = 0.5, color = "darkblue", fill = "lightblue") +
  labs(y = "Histogram", x = "Residuals") +
  theme(text = element_text(size=15))
pResidNorm
#' Looking more normal than model 1, but still not as good as it could be.

#' Based on the scaled quantile residuals, we can probably assume that the plots
#' of residuals vs covariates will have issues, but let's just see one or two.

#' Plot residuals vs. age.years
pResidvsAge <- ggplot(data = schisto.pos.sf, aes(y = resid.mod2, x = age.years)) + 
  geom_point(shape = 1, size = 1) +
  geom_smooth() +
  labs(x = "Age (years)", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsAge
#' Not ideal

#' Plot residuals vs. EVI
pResidvsEVI <- ggplot(data = schisto.pos.sf, aes(y = resid.mod2, x = jan_apr_23_evi.std)) + 
  geom_point(shape = 1, size = 1) +
  geom_smooth() +
  labs(x = "EVI (std)", y = "Residuals") +
  theme(text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2)
pResidvsEVI
#' Not too bad, but not great either.

#' Step 4: Check for spatial dependency
#' There are 3 ways we can do this:
#'  1. Plot residuals vs spatial locations and look for patterns
#'  2. Apply a Moran's I test on the residuals
#'  3. Make a variogram of the residuals
#' We'll do options 1 and 3 here.

#' Option 1: Plot residuals vs spatial locations

#' Define the color and point size for the plot
#' Remember that scaled quantile residuals are centered on 0.5, while Pearson 
#' residuals are centered around 0.
schisto.pos.sf$MyCol <- ifelse(schisto.pos.sf$resid.mod2 >= 0, "red", "blue")
schisto.pos.sf$MySize <- rescale(abs(schisto.pos.sf$resid.mod2), to = c(0,10))

#' Plot it
ggplot(data = kruger.km) +
  geom_sf(fill = "transparent") +
  geom_sf(data = schisto.pos.sf,
          aes(col = MyCol,
              size = MySize),
          alpha = 0.3) +
  scale_color_identity() +
  scale_size_continuous(range = c(1,3)) +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Residuals")
#' Hmm, hard to tell.

#' Option 3: Make a variogram
#' First, make a dataframe with the residuals and capture locations
mod2data <- data.frame(E1 = schisto.pos.sf$resid.mod2,
                       Xkm = schisto.pos.sf$Xkm,
                       Ykm = schisto.pos.sf$Ykm)

#' Convert to an `sf` object
mod2data_sf <- st_as_sf(x = mod2data,
                        coords = c("Xkm", "Ykm"),
                        crs = NA) #' CRS doesn't matter for this, so no need to specify

#' Call the `variogram()` function from the package `gstat`
V2 <- variogram(E1 ~1,
                data = mod2data_sf,
                cutoff = 100,
                cressie = TRUE)

#' And plot
ggplot(data = V2, aes(x = dist, y = gamma)) +
  geom_point() +
  geom_smooth(se = FALSE, span = 1) +
  labs(x = "Distance (km)", y = "Semi-variogram") +
  theme(text = element_text(size = 15))
#' This looks a lot better than without the SRF.


# Section 9: Model Comparison

#' We concluded that adding in the SRF at least removed some spatial dependency
#' in the residuals, but let's compare DICs and WAICs to check that the SRF
#' model is indeed better. Keep in mind that based off of the residuals, even
#' the SRF model needs serious work.

DICs  <- c(mod1$dic$dic, mod2$dic$dic)
WAICs <- c(mod1$waic$waic, mod2$waic$waic)

Results <- data.frame(Models = c("lm",
                                 "lm + SRF"),
                      DIC = DICs,
                      WAIC = WAICs)
colnames(Results) <- c("  Models", 
                       "       DIC", 
                       "       WAIC")
Results

#' Based off of this, the model without the SRF has a higher DIC but a lower WAIC,
#' so it's difficult to conclude anything. 

#' It's probably not worth it, but let's plot the SRF and spatial dependency from
#' model 2.

#' Posterior distribution of the range.
Range     <- spde.posterior(mod2, "SRF", what = "range")
plot(Range)


#' This is the imposed spatial dependence.
MaternCor <- spde.posterior(mod2, "SRF", what = "matern.correlation")
pB <- plot(MaternCor) +
  labs(x = "Distance (km)", y = "Correlation") +
  coord_cartesian(ylim = c(0, 1))
pB

#' Define strong correlation as correlation between 1 - 0.8.
#' Define moderate correlation as correlation between 0.8 - 0.5.
#' Define weak correlation as correlation between 0.5 - 0.1.
#' Define diminishing correlation as correlation smaller than 0.1.


#' Would you like to add the 0.8, 0.5 and 0.1 thresholds?
#' The code is slightly intimidating, feel free to just run it.
gb   <- ggplot_build(pB)
dfB  <- gb$data[[1]][, c("x","y")]
dfB  <- dfB[order(dfB$x), ]
dfB2 <- dfB[order(dfB$y), ]  # y ascending, x matched


#' Helper: distance at which correlation hits a given level (linear interpolation)
x_at <- function(level) {
  approx(x = dfB2$y, y = dfB2$x, xout = level, rule = 2)$y
}

thr   <- c(0.8, 0.5, 0.1)
xthr  <- sapply(thr, x_at)
marks <- data.frame(level = thr, x = xthr)

#' Make the guide segments + labels
xmax <- max(dfB$x)
p.Matern <- pB +
  geom_segment(data = marks, 
               aes(x = 0, xend = x, y = level, yend = level), 
               linetype = "dotted") +
  geom_segment(data = marks, 
               aes(x = x,  xend = x, y = 0, yend = level), 
               linetype = "dotted") +
  geom_point(data = marks, 
             aes(x = x, y = level), 
             size = 1.2) +
  annotate("text", 
           x = xthr, y = 0, 
           label = paste0("~", round(xthr, 1)),
           vjust = 1.4, 
           size = 3) +
  annotate("text", 
           x = 0.97 * xmax, y = 0.90, 
           hjust = 1, 
           size = 3.2,
           label = "Strong: 1–0.8") +
  annotate("text", 
           x = 0.97 * xmax, 
           y = 0.65, 
           hjust = 1, 
           size = 3.2,
           label = "Moderate: 0.8–0.5") +
  annotate("text", 
           x = 0.97 * xmax, 
           y = 0.30, 
           hjust = 1, 
           size = 3.2,
           label = "Weak: 0.5–0.1") +
  annotate("text", 
           x = 0.97 * xmax, 
           y = 0.06, 
           hjust = 1, 
           size = 3.2,
           label = "Diminishing: < 0.1") +
  theme_minimal() +
  theme(text = element_text(size = 15)) 
p.Matern
#' This says that the spatial autocorrelation drops below 0.8 at ~5.3km, below
#' 0.5 at ~11.8km, and below 0.1 at ~30.2km.

#' Now let's plot the SRF
InnerHull_sfc <- fm_as_sfc(inner.hull) 
st_crs(InnerHull_sfc) <- st_crs(schisto.pos.sf) 

#' We'll use the outer hull, so that we have all of KNP covered
OuterHull_sfc <- st_geometry(outer_sf_simple)

NewData <- fm_pixels(mesh1, 
                     dims = c(250, 250),
                     mask = OuterHull_sfc)
Pred <- predict(object = mod2, 
                newdata = NewData, 
                formula = ~ SRF)

#' Some nice colors from the inlabru website:
#' https://inlabru-org.github.io/inlabru/articles/random_fields.html.
ColSc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE))}

#' Plot the spatial random field.
ggplot() +
  gg(Pred, geom = "tile") +
  geom_sf(data = kruger.km,
          fill = "transparent", 
          col = "darkred") +
  ggtitle("Posterior mean SRF") +
  labs(x = "Easting (km)", y = "Northing (km)") + 
  ColSc(Pred$mean)

#' Recall that the output of the figure is information in CAA_i that **cannot**
#' be explained by the covariates and is spatial in nature!

#' Refer to the `mod2_exercise4_inlabru_lm_spatial.R` file from the Highland
#' Statistics inla course for code on visualizing the u_i values for the SRF.