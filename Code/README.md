# Code

This folder contains all the code for the NBDA_SBdolphins repository. The following describes the analysis steps in the `NBDA_Bayes_sd.R` file. 

## Data analysis process
<img src="https://github.com/user-attachments/assets/910b3c00-0a76-4735-b39b-6cb383488591" align="middle" width="500px"/>

## Supplemental Step 1: Data Wrangling (code available upon request)

I start by cleaning up and combining data from 1995-2014. I then evaluating the data using plots and diagnosing code.

## PART 1: Create Networks

I create the horizontal, vertical and ecological networks for each inter-event period. 

## PART 2: Create acquisition data for model input

Here I create all aquisition data.

## PART 3: Rund the models

I run a multi-network-based diffusion analysis using a Markov chain Monte Carlo (MCMC) sampler under a Bayesian statistical framework. 

Process model:
$SRI_{i,j,p} ~ (u_i+u_j)+β1_{Sex_{i,j}}+β2_{Age_{i,j}}+β3_{HRO_{i,j,p}}+β4_{HC_{i,j,p}}+β5_{During HAB_{i,j,p}}+β6_{After HAB_{i,j,p}}+β7_{HC_{i,j,p}}*During HAB+β8_{HC_{i,j,p}}*After HAB$

Observation model:
$SRI_{i,j,p} ~ Beta[μ_{SRI_{i,j,p}},ϕ_{SRI_{i,j,p}}]$

## PART 4: Summary Outputs

I create summary outputs for the model.
