# Estimating the effect of social inequalities on the mitigation of COVID-19 across communities in Santiago de Chile
Data and code for the paper "Estimating the effect of social inequalities in the mitigation of COVID-19 across communities in Santiago de Chile"

![alt text](https://github.com/ngozzi/covid19-santiago/blob/master/santiago.png)

https://www.medrxiv.org/content/10.1101/2020.10.08.20204750v1

## Abstract
We study the spatio-temporal spread of SARS-CoV-2 in Santiago de Chile using anonymized mobile phone data from 1.4 million users, 22% of the whole population in the area, characterizing the effects of non-pharmaceutical interventions (NPIs) on the epidemic dynamics. We integrate these data into a mechanistic epidemic model calibrated on surveillance data. As of August 1, 2020, we estimate a detection rate of 102 cases per 1,000 infections (90% CI: [95 - 112 per 1,000]). We show that the introduction of a full lockdown on May 15, 2020, while causing a modest additional decrease in mobility and contacts with respect to previous NPIs, was decisive in bringing the epidemic under control, highlighting the importance of a timely governmental response to COVID-19 outbreaks. We find that the impact of NPIs on individuals' mobility correlates with the Human Development Index of comunas in the city. Indeed, more developed and wealthier areas became more isolated after government interventions and experienced a significantly lower burden of the pandemic. The heterogeneity of COVID-19 impact raises important issues in the implementation of NPIs and highlights the challenges that communities affected by systemic health and social inequalities face adapting their behaviors during an epidemic.

## Data 
The data folder contains the commuting matrices across comunas in the three different regimes (baseline, soft and full lockdown) and the contacts reduction parameters in different comunas during the soft and full lockdown.

## Code
### Compile and run
The code folder contains the C++ file needed to run the epidemic model. To compile and generate an executable file, open the code directoty in terminal and type: 
g++ -std=c++11 main.cpp include/sampler.h include/sampler.cpp include/Parser.h 
This will generate an executable file called "a.out". To run it, just type "./a.out" in ther terminal window. 

### Input
The model first parses the input files which are in the "input" folder. These files are in json format and they contain the commuting matrices (baseline, soft and full lockdown), the population files, the contacts matrices and the contacts reduction parameters. The parsing is done with the following library: https://github.com/nlohmann/json. To change the epidemiological parameters (i.e. R0), modify directly the main.cpp file (before running again, remember to recompile).

### Output
The model generates a "results.txt" file in the output directory. This file contains the evolution in time of the number of S, E, I, R individuals in different comunas and different age groups. For example the column "I_10_5", contains the the time series of the number of infected for comuna 10 in age group 5. 

## Other Models
This repository includes also the code for other models. In particular, in the folder "models" there are two additional models:
- SLPIAR.cpp: this model extends the metapopulation one presented in the file main.cpp including a more complex disease dynamics, with asymptomatic and presymptomatic infectious individuals.
- no_mobility.cpp: this model represent the comunas as disconnected subpopulations. The unfolding of the spreading is thus represented here without mobility coupling, but only with contacts reduction parameters estimated via mobile phone data.

The folder "single-pop" includes the python code for two models in which we disregard the organization of Santiago in connected subpopulations (comunas) and we represent it as a simple single population. The two models differ for the kind of data used to estimates contact reduction. In the first one, we use our mobile phone data to infer the overall contacts reduction in Santiago following the establishment of restrictions. In the second one, we use instead data from the [COVID-19 Google Community Mobility Report](https://www.google.com/covid19/mobility/) and from the [Oxford COVID-19 Government Response Tracker](https://www.bsg.ox.ac.uk/research/research-projects/covid-19-government-response-tracker).

