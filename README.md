# Schisto_Analysis
Spatial patterns of schistosoma infection in African buffalo in KNP using INLA-SPDE

The purpose of the `test_schisto_spde.Rmd` file is to see if we can resolve our issue of lack of spatial variability by only including 1 buffalo per unique capture location. For the purposes of this analysis, I only included one environmental variable (elevation).
The `schisto_bym_1.Rmd` file was a preliminary analysis using INLA-BYM and BYM2 models to see if there were any *sections* of KNP that would be predicted to have a higher likelihood of infection, but these analyses do *not* include any environmental covariates.

The current issues are:
1. Trying to run the spatial model with a reasonable resolution (mesh triangles, in this case), but without overpowering the model. I think right now there are ~190k nodes for only 64 observations, which I understand is quite a lot.
2. Getting the mesh and the boundary to project in the right space - for some reason I cannot get them both to project in the correct hemisphere.

Of course, if you see anything else in the code that I've done incorrectly, please feel free to edit! I included the model code in the .Rmd file, rather than running them in the background, just so I could troubleshoot a bit better.

Thank you Henri!