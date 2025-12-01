R sciprts for paper: Bolam, S.G., Cooper, K.M., and 
Downie, A-L. Developing an Ecological Risk-Based Approach to Facilitate
Licensing Offshore Wind Development. Ecosphere.

In this study we create a Combined Risk layer based on 3 element layers for: 
Biodiversity, Sensitivity and Assemblage Rarity.

PART 1: SENSITIVITY (VERSION 1.0. 2025; see sensitivity.R) is used to develop the  
Sensitivity layer. Whilst the code includes lines for generating a random forest
Sensitivity model, this is intended only as a quick look see. The final 
Sensitivity model, and associated confidence layer should be created using code
from Risk_ContinuousVariablesModel_2025.R, using input data (i.e. point sample 
sensitivity scores) generated in this file.

PART 2: ASSEMBLAGES (VERSION 3.0. 2025; see assemblages.R)' is used to develop the  
Assemblages layer. Whilst the code includes lines for generating a random forest
Assemblages model, this is intended only as a quick look see. The final Assemblages
model, and associated confidence layer should be created using code from 
Risk_ClusterModel_2025.R, with input data (i.e. point sample assemblages) generated
in this file.

The Biodiversity layer was created in a sperate study, with relevant .R files found
here: https://github.com/keithmcooper/OneBenthicBiodiversity

PART 3: RISK (risk.R) is used to develop the final Combined Risk layer,
together with an associated confidence layer. Input data are raster files for 
Sensitivity (see PART 1: SENSITIVITY;sensitivity.R), Biodiversity 
(see https://github.com/keithmcooper/OneBenthicBiodiversity), and Assemblage Rarity
(see PART 2: ASSEMBLAGES; assemblages.R).
