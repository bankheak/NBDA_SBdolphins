# Code

This folder contains all the code for the NBDA_SBdolphins repository. The following describes the analysis steps in the `NBDA_Bayes_sd.R` file. 

## Data analysis process
<img width="1461" height="622" alt="NBDA_workflow" src="https://github.com/user-attachments/assets/7312c167-ed93-474a-b76c-74a8117aadfa" />

## Supplemental Step 1: Data Wrangling (code available upon request)

I start by cleaning up and combining data from 1995-2014. I then evaluating the data using plots and diagnosing code.

## PART 1: Create Networks

I create the horizontal, vertical and ecological networks for each inter-event period. 

Within the horizontal model code I add a column for inter-event period that is calculated based if new individuals aquired the target behavior for each month. 
The group-by-individual data frame is then split into a list based on each aquisition event and an association matrix is calculated for each element in the list.
Observations of individuals after they aquired the target behavior were dropped from the observation data used to calculate association.

Dyads in the vertical network were simply given a 1 if they were a mother-calf pair and a 0 if they were not.

The ecological network was created using the kernel density home range overlap estimates of individuals. These were split by year in order to get a more robust reading
of individuals' home ranges.

## PART 2: Create acquisition data for model input

Here I create all aquisition event data, edge lists, individual-level variables, and weight estimates to input as predictors and response variables for the NBDA model.
I started by making sure that the networks were all the same dimensions based on the number of inter-event periods and included the same individuals. I then combined these
networks into one edgelist data frame that included time and network estimates for each dyad. 

To create the event data I calculated each individual's aquisition time, 
the overall end time based on the last aquisition event. Individuals that never aquired the target behavior were assigned the end time + 1.

I then created the constant individual-level variable (ILV_c) and varying individual-level variable (ILV_v) dataframes. These included sex as constant
and age as varying.

Finally, the proportion of time each individual spent engaging in human-centric behavior was calculated and added as a static weight in the data list.

## PART 3: Create acquisition data for model input

I read in the event_data, edge_list, ILV_c, ILV_tv, and HI_matrix wrangled from the previous step and combined it into a import_user_STb() data list to be run in the model.

## PART 4: Run the model

I ran a multi-network-based diffusion analysis using a Markov chain Monte Carlo (MCMC) sampler under a Bayesian statistical framework. I used a test model to run the raw data and found that 
the model detected N_veff = 0, which meant that the likihood provided no information, this lead me to change the baseline learning rate to be positive and bounded and change the priors.
I did this by working in the STAN.model code which were saved for each behavior. This was then ran in fit_STb().


NBDA model:

$\lambda_i(t) \sim \lambda_0(t)(1-z_i(t))\left(e^{\Gamma_i(t)}\sum_n S_n\sum_{j=1}^{N}w_j(t)a_{n,ij}(t)z_j(t)+e^{\beta_i(t)}\right)$

$\beta_i \sim \sum_{k=1}^{V} \beta_k x_{k,i}$

$\Gamma_i \sim \sum_{k=1}^{V} \gamma_k x_{k,i}$

where λi(t) is the rate at which individual i acquires the target behavior as a function of time, λo(t) is a baseline rate function, zi(t) is the ‘status’ of individual i at time t (1 = informed; 0 = naïve), N is the number of individuals in the population, w_j is the transmission weight of the rate at which each individual performed human-centric foraging behavior during difusion, n is the number of networks and aij indicates the connection strength from j to i from the social networks, xk,i is the value of the kth variable for individual i, βk is the coefficient of the effect of variable k (sex, age and HAB exposure) on asocial learning, and γk is the coefficient of the effect of variable k on social transmission. The key model output will be the relative strength of social transmission, s, the value of which is estimated when the model is fitted to the data.

## PART 5: Summary Outputs

I first ran posterior predictive checks to make sure the observed data lines up with model predictions.

I then plotted model predictor effect sizes separated by social and individual ILV versus network effects.

Finally I extracted hazard rates for target behaviors. To do this I first extracted posterior distributions and the social and individual inputs for each network over each inter-event period. I then used these estimates to compute hazard rates per period, across all posterior draws. I then took the summary of these and plotted them to understand the difference in social versus individual learning over time. Finally, I found the trial when social learning overtook individual learning and related it to the study year.
