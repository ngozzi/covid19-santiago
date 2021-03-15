# libraries
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from datetime import date as dt_date 
import uuid
import pickle as pkl 
import gzip
from numba import jit
import warnings
warnings.filterwarnings("ignore")


# n. of compartments and age groups
ncomp = 4
nage  = 16


# epidemiological parameters
IFR = [0.00161 / 100, # 0-4  
       0.00161 / 100, # 5-9
       0.00695 / 100, # 10-14
       0.00695 / 100, # 15-19 
       0.0309  / 100, # 20-24 
       0.0309  / 100, # 25-29
       0.0844  / 100, # 30-34
       0.0844  / 100, # 35-39
       0.161   / 100, # 40-44 
       0.161   / 100, # 45-49 
       0.595   / 100, # 50-54 
       0.595   / 100, # 55-59 
       1.93    / 100, # 60-64
       1.93    / 100, # 65-69 
       4.28    / 100, # 70-74 
       6.04    / 100] # 75+


mu  = 1 / 2.5
eps = 1 / 4.0


# fit dates
start_date = datetime(2020, 3, 1)
end_date   = datetime(2020, 8, 2)
date       = start_date




def import_basin():
    
    basin_dict = dict()
    
    # import restriction data
    basin_dict["other_loc_restr"] = pd.read_csv("../data/restrictions/other_loc.csv")
    basin_dict["work_restr"]      = pd.read_csv("../data/restrictions/work.csv")
    basin_dict["school_restr"]    = pd.read_csv("../data/restrictions/school.csv")
    basin_dict["home_restr"]    = pd.read_csv("../data/restrictions/home.csv")
    basin_dict["overall_restr"]    = pd.read_csv("../data/restrictions/overall.csv")
    
    # import contacts matrices
    basin_dict["other_loc_contacts"] = pd.read_csv("../data/contacts-matrix/chile_other_locations.csv", header=None).values
    basin_dict["work_contacts"]      = pd.read_csv("../data/contacts-matrix/chile_work.csv", header=None).values
    basin_dict["school_contacts"]    = pd.read_csv("../data/contacts-matrix/chile_school.csv", header=None).values
    basin_dict["home_contacts"]      = pd.read_csv("../data/contacts-matrix/chile_home.csv", header=None).values
    
    # import initial conditions
    with open("../data/initial-conditions/exposed.pkl", "rb") as file:
        basin_dict["exposed_ic"] = pkl.load(file)
    with open("../data/initial-conditions/infectious.pkl", "rb") as file:
        basin_dict["infectious_ic"] = pkl.load(file)
        
    # import Nk
    basin_dict["Nk"] = pd.read_csv("../data/demographic/pop_5years.csv")["total"].values
    
    # import epi data
    basin_dict["epi_data"] = pd.read_csv("../data/epidemiological-data/epi_data.csv")
    
    return basin_dict


def update_contacts_TEF(basin_dict, date):

    # get baseline contacts matrices
    home      = basin_dict["home_contacts"]
    work      = basin_dict["work_contacts"]
    school    = basin_dict["school_contacts"]
    other_loc = basin_dict["other_loc_contacts"]
        
    # get overall reductions   
    r1 = 0.5621746959575077
    r2 = 0.4516456264414326

    if date <= datetime(2020, 3, 16):
        omega = 1
    elif date > datetime(2020, 3, 16) and date < datetime(2020, 5, 15):
        omega = r1
    else:
        omega = r2
    # contacts matrix with reductions
    C = omega * (home +  school + work + other_loc)
    
    return C


def update_contacts_GOOGLE(basin_dict, date):

    # get baseline contacts matrices
    home      = basin_dict["home_contacts"]
    work      = basin_dict["work_contacts"]
    school    = basin_dict["school_contacts"]
    other_loc = basin_dict["other_loc_contacts"]
        
    # get year-week
    if date.isocalendar()[1] < 10:
        year_week = str(date.isocalendar()[0]) + "-0" + str(date.isocalendar()[1])
    else:
        year_week = str(date.isocalendar()[0]) + "-" + str(date.isocalendar()[1])
        
    # get work / other_loc / school reductions   
    work_reductions           = basin_dict["work_restr"]
    other_loc_reductions      = basin_dict["other_loc_restr"]
    home_reductions           = basin_dict["home_restr"]
    school_reductions         = basin_dict["school_restr"]
    school_reductions["date"] = pd.to_datetime(school_reductions["date"])

    omega_h = home_reductions.loc[home_reductions.year_week==year_week]["home_red"].values[0]
    omega_w = work_reductions.loc[work_reductions.year_week==year_week]["work_red"].values[0]
    omega_o = other_loc_reductions.loc[other_loc_reductions.year_week==year_week]["oth_red"].values[0]
    omega_s = (3 - school_reductions.loc[school_reductions.date==date]["C1_School closing"].values[0]) / 3
    
    overall_reductions = basin_dict["overall_restr"]
    overall_reductions["date"] = pd.to_datetime(overall_reductions["date"])
    omega = overall_reductions.loc[overall_reductions["date"]==date].values[0][1]
   
    # contacts matrix with reductions
    C = (omega_h * home) + (omega_s * school) + (omega_w * work) + (omega_o * other_loc)

    return C



@jit(fastmath=True)
def stochastic_SEIRD(Cs, basin_dict, initial_conditions, beta, dates):

    """
    This function simulates a stochastic SEIR model
    """

    # initial conditions
    T = len(dates)
    compartments = np.zeros((ncomp, nage, T))
    compartments[:,:,0] = initial_conditions

    # simulate 
    for i in range(T - 1): 

        C = Cs[dates[i]]

        # next step solution 
        next_step = np.zeros((ncomp, nage))

        # iterate over ages
        for age1 in range(nage):

            # compute force of infection
            # S:0 L:1 I:2 R:3
            force_inf = np.sum(beta * C[age1, :] * compartments[2, :, i] / basin_dict["Nk"])

            # S -> L
            if force_inf == 0:
                new_latent = 0
            else:
                new_latent  = np.random.binomial(compartments[0, age1, i], force_inf)

            # L -> I 
            new_infected  = np.random.binomial(compartments[1, age1, i], eps)

            # I -> R
            new_recovered = np.random.binomial(compartments[2, age1, i], mu)

            # update next step solution
            next_step[0, age1] = compartments[0, age1, i] - new_latent                       # S
            next_step[1, age1] = compartments[1, age1, i] + new_latent   - new_infected      # L
            next_step[2, age1] = compartments[2, age1, i] + new_infected - new_recovered     # I
            next_step[3, age1] = compartments[3, age1, i] + new_recovered                    # R

        # update solution at the next step 
        compartments[:,:,i+1] = next_step

    return compartments



def simulate(Cs, basin_dict, initial_conditions, R0, Delta, dates):

    """
    This function runs the SEIR model
    """

    # get beta from Rt w.r.t to initial_date
    beta = get_beta(R0, mu, Cs[dates[0]])

    # simulate
    compartments = stochastic_SEIRD(Cs, basin_dict, initial_conditions, beta, dates)

    # compute deaths
    deaths = compute_deaths(compartments, Delta)

    return compartments, deaths



@jit(fastmath=True)
def compute_deaths(compartments, Delta):

    """
    This function computes the number of daily deaths
    """

    # initialize deaths
    nage, T = compartments.shape[1], compartments.shape[2]
    deaths  = np.zeros((nage, T))

    # get recovered of this step
    recovered = np.diff(compartments[3], axis=1)

    for age in range(nage):
        deaths[age] = np.concatenate(([0], np.random.binomial(recovered[age].astype(int), IFR[age])))
        deaths[age] = np.roll(deaths[age], shift=Delta)
        deaths[age][0:Delta] = 0

    return deaths


def get_beta(R0, mu, C):

    """
    This functions return beta for a SEIR model with age structure
    """

    # get seasonality adjustment
    return R0 * mu / (np.max([e.real for e in np.linalg.eig(C)[0]]))



# import data on the basin
basin_dict = import_basin()
    
    
# real deaths
tot_deaths = basin_dict["epi_data"]
tot_deaths.index = pd.to_datetime(tot_deaths["fecha"])
tot_deaths = tot_deaths.resample("W").sum().loc[tot_deaths.resample("W").sum().index <='2020-08-02']["count"].values
    
    
# initial conditions
initial_conditions = np.zeros((ncomp, nage))
age_cols = ['0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', 
            '60-64', '65-69', '70-74', '75+']
for age_id, age_str in zip(range(16), age_cols):
    # (S: 0, L1: 1, L2: 2, I1: 3, I2: 4, R1: 5, R2: 6)
    initial_conditions[1, age_id] = int(basin_dict["exposed_ic"][age_str]) 
    initial_conditions[2, age_id] = int(basin_dict["infectious_ic"][age_str]) 
    initial_conditions[3, age_id] = int(0)
    initial_conditions[0, age_id] = int(basin_dict["Nk"][age_id] - initial_conditions[1, age_id] - initial_conditions[2, age_id])



# pre-compute contacts matrices over time
Cs_TEF, Cs_GOOGLE = {}, {}
date, dates = start_date, [start_date]
for i in range((end_date - start_date).days - 1):
    Cs_TEF[date] = update_contacts_TEF(basin_dict, date)
    Cs_GOOGLE[date] = update_contacts_TEF(basin_dict, date)
    date += timedelta(days=1)
    dates.append(date)


# sample run
R0 = 2.66
Delta = 17

# TEF
compartments_TEF, deaths_TEF = simulate(Cs_TEF, basin_dict, initial_conditions, R0, Delta, dates)

# GOOGLE / OXFORD
compartments_GOOGLE, deaths_GOOGLE = simulate(Cs_GOOGLE, basin_dict, initial_conditions, R0, Delta, dates)
