# SPSP 0.2.0

## New features

1. Add the supports to the SCAD-penalized and MCP-penalized regression model through the `fitfun.SP = SCAD.ncvreg`, `fitfun.SP = MCP.ncvreg`. The previous established methods are `ridge.glmnet`. `lasso.glmnet`,
`adalasso.glmnet`, and `adalassoCV.glmnet`.

# SPSP 0.1.1

## patch changes

1. Fix the function corresponding to the `intercept `argument. The intercept argument controls whether to include the intercept, but since the SPSP algorithm estimates the intercept separately. It needs to be dealt with separately.

2. The `standardize` only control the standardization of covariates.

# SPSP 0.1.0

## Major changes

1. Implement the Selection by Partitioning the Solution Paths procedure in the SPSP() and SPSP_step() functions;

2. Add fitfun.SP functions to obtain the solution paths for the SPSP algorithm;

3. Include one high-dimensional data for illustration purpose.



