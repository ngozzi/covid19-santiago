#include <iostream>
#include <vector>
#include <cstdlib>
#include "./include/sampler.h"
#include "./include/Parser.h"

using namespace std;


/**
 * Compute beta given R0
 */
double get_beta(double R0, double alpha, double gamma, double frac_asympt, double mu, double eigenV)
{
    double beta = R0 / (eigenV * (alpha / gamma + (1 - frac_asympt) / mu + alpha * frac_asympt / mu));
    return beta;
}



/**
 * Multiply matrix by scalar
 */
vector<vector<double>> scalarProductMat(vector<vector<double>> C, double r, int K)
{
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            C[i][j] = C[i][j] * r;
    return C;
}


/**
 * Sum two matrix with same dimension
 */
vector<vector<double>> sumMat(vector<vector<double>> C1, vector<vector<double>> C2, int K)
{
    vector<vector<double>> C(K, vector<double>(K, 0.0));
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j ++)
            C[i][j] = C1[i][j] + C2[i][j];
    return C;
}


/**
 * Get sigma_i for population i
 */
double get_sigma(int j, vector<vector<double>> &sigmas)
{
    double sigma = 0.0;
    for (int i = 0; i < sigmas[j].size(); i++)
        sigma += sigmas[j][i];
    return sigma;
}


/**
 * Get effective N for a specific age group
 */
double get_Nk_eff(int pop_idx, int age_idx, double tau, vector<vector<double>> &Nk, vector<vector<double>> &sigmas, vector<double> &sigmas_j)
{
    double Neff = Nk[pop_idx][age_idx] / (1 + sigmas_j[pop_idx] / tau);
    for (int i = 0; i < Nk.size(); i++)
        Neff += Nk[i][age_idx] / (1 + sigmas_j[i] / tau) * sigmas[i][pop_idx] / tau;
    return Neff;
}


/**
 * Get single lambda
 */
double get_lambda(int pop_idx, int age_idx, double tau, double beta, double alpha, vector<vector<double>> &Nk_eff, vector<vector<double>> &P, vector<vector<double>> &I, vector<vector<double>> &A, vector<vector<double>> &C, vector<vector<double>> &sigmas, vector<double> &sigmas_j)
{
    double lambda = 0.0;
    for (int i = 0; i < I.size(); i++)
        if (i == pop_idx)
            for (int k = 0; k < I[0].size(); k++)x
                lambda += (C[age_idx][k] / Nk_eff[pop_idx][k]) * ((I[pop_idx][k] + alpha * (P[pop_idx][k] + A[pop_idx][k])) / (1 + sigmas_j[pop_idx] / tau));
        else
            for (int k = 0; k < I[0].size(); k++)
                lambda += (C[age_idx][k] / Nk_eff[pop_idx][k]) * ((I[i][k] + alpha * (P[i][k] + A[i][k])) / (1 + sigmas_j[i] / tau)) * (sigmas[i][pop_idx] / tau);
    lambda = lambda * beta;
    return lambda;
}


/**
 * Get total lambda
 */
double get_lambda_tot(int pop_idx, int age_idx, double tau, double beta, vector<vector<double>> &Nk_eff, vector<vector<double>> &I, vector<vector<vector<double>>> &C, vector<vector<double>> &sigmas, vector<double> &sigmas_j)
{
    double lambda_tot = 0.0;
    double lambda_ji = 0.0;
    double lambda_jj = get_lambda(pop_idx, age_idx, tau, beta, Nk_eff, I, C[pop_idx], sigmas, sigmas_j);
    lambda_tot += lambda_jj / (1 + sigmas_j[pop_idx] / tau);

    for (int i = 0; i < I.size(); i++)
    {
        lambda_ji = get_lambda(i, age_idx, tau, beta, Nk_eff, I, C[i], sigmas, sigmas_j);
        lambda_tot += lambda_ji / (1 + sigmas_j[pop_idx] / tau) * (sigmas[pop_idx][i] / tau);
    }
    return lambda_tot;
}


int main(int argc, char *argv[])
{

    // number of comunas, age groups, and simulations per parameters
    int Npop = 39;
    int K = 16;
    int Nsim = 200;

    // parameters
    double mu = 1 / 2.5;
    double eps = 1 / 3.7;
    double gamma = 1 / 1.5;
    double frac_asympt = 0.4;
    double alpha = 0.55;
    double tau = 3.0;

    vector<double> betas;
    for (double R0 = 2.0; R0 < 4.02; R0 = R0 + 0.02)
        betas.push_back(get_beta(R0, alpha, gamma, frac_asympt, mu, 16.204308331681283));

    for (int b = atoi(argv[1]); b < atoi(argv[2]); b++)
    {
        double beta = betas[b];
        for (int l = 0; l < Nsim; l++)
        {
            // parser
            Parser parser = Parser("/mnt/beegfs/home/ng9035y/chile_rebuttal/SLPIAR/input/config" + to_string(l) + ".json");

            // commuting
            vector<vector<double>> sigmas = parser.parse_commuting();
            vector<double> sigmas_j;
            for (int i = 0; i < Npop; i++)
                sigmas_j.push_back(get_sigma(i, sigmas));

            // compartments
            vector<vector<double>> S = parser.parse_compartments("S");
            vector<vector<double>> L = parser.parse_compartments("L");
            vector<vector<double>> P = parser.parse_compartments("P");
            vector<vector<double>> I = parser.parse_compartments("I");
            vector<vector<double>> A = parser.parse_compartments("A");
            vector<vector<double>> RI = parser.parse_compartments("RI");
            vector<vector<double>> RA = parser.parse_compartments("RA");
            vector<vector<double>> Nk = parser.parse_compartments("N");


            // contacts (home, other)
            vector<vector<double>> C1 = parser.parse_contacts(1);
            vector<vector<double>> C2 = parser.parse_contacts(2);
            vector<vector<vector<double>>> C;
            for (int i = 0; i < Npop; i++)
                C.push_back(sumMat(C1, C2, K));

            vector<double> r1 = parser.parse_r(1);
            vector<double> r2 = parser.parse_r(2);

            // compute N effective k
            vector<vector<double>> Nk_eff(Npop, vector<double>(K, 0.0));
            for (int i = 0; i < Npop; i++)
                for (int k = 0; k < K; k++)
                    Nk_eff[i][k] = get_Nk_eff(i, k, tau, Nk, sigmas, sigmas_j);

            // next step data
            vector<vector<double>> Snext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> Lnext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> Pnext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> Inext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> Anext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> RInext(Npop, vector<double>(K, 0.0));
            vector<vector<double>> RAnext(Npop, vector<double>(K, 0.0));
            double newL = 0.0;
            double newP = 0.0;
            double newIA = 0.0;
            double newI = 0.0;
            double newA = 0.0;
            double newRI = 0.0;
            double newRA = 0.0;
            double lambda = 0.0;

            // write results header
            ofstream resFile("/mnt/beegfs/home/ng9035y/chile_rebuttal/SLPIAR/output/results_" + to_string(b) + "_" + to_string(l) + ".txt");
            for (int i = 0; i < Npop; i++)
                for (int k = 0; k < K; k++)
                    resFile << "P_" << to_string(i) << "_" << to_string(k) << "," << "I_" << to_string(i) << "_" << to_string(k) << "," << "A_" << to_string(i) << "_" << to_string(k) << ","<< "RI_" << to_string(i) << "_" << to_string(k) << "," << "RA_" << to_string(i) << "_" << to_string(k) << "," ;
            resFile << "\n";

            // simulate
            for (int t = 60; t < 250; t++)
            {
                if (t == 75) // restrictions
                {
                    // import new commuting, recompute C and Nk_eff
                    Parser parser = Parser("/mnt/beegfs/home/ng9035y/chile_rebuttal/SLPIAR/input/config_restr.json");

                    C.clear();
                    for (int i = 0; i < Npop; i++)
                        C.push_back(sumMat(scalarProductMat(C1, r1[i], K), scalarProductMat(C2, r1[i], K), K));

                    sigmas.clear();
                    sigmas_j.clear();
                    sigmas = parser.parse_commuting();
                    for (int i = 0; i < Npop; i++)
                        sigmas_j.push_back(get_sigma(i, sigmas));

                    for (int i = 0; i < Npop; i++)
                        for (int k = 0; k < K; k++)
                            Nk_eff[i][k] = get_Nk_eff(i, k, tau, Nk, sigmas, sigmas_j);
                }

                else if (t == 136)  // lockdown
                {
                    // import new commuting, recompute C and Nk_eff
                    Parser parser = Parser("/mnt/beegfs/home/ng9035y/chile_rebuttal/SLPIAR/input/config_restr1.json");

                    C.clear();
                    for (int i = 0; i < Npop; i++)
                        C.push_back(sumMat(scalarProductMat(C1, r2[i], K), scalarProductMat(C2, r2[i], K), K));

                    sigmas.clear();
                    sigmas_j.clear();
                    sigmas = parser.parse_commuting();
                    for (int i = 0; i < Npop; i++)
                        sigmas_j.push_back(get_sigma(i, sigmas));

                    for (int i = 0; i < Npop; i++)
                        for (int k = 0; k < K; k++)
                            Nk_eff[i][k] = get_Nk_eff(i, k, tau, Nk, sigmas, sigmas_j);
                }


                for (int i = 0; i < Npop; i++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        // S -> L
                        lambda = get_lambda_tot(i, k, tau, beta, alpha, Nk_eff, P, I, A, C, sigmas, sigmas_j);
                        newL = binomial(S[i][k], lambda);

                        // L -> P
                        newP = binomial(L[i][k], eps);

                        // P -> I/A
                        newIA = binomial(P[i][k], gamma);
                        newA  = binomial(newIA, frac_asympt);
                        newI  = newIA - newA

                        // I -> RI
                        newRI = binomial(I[i][k], mu);

                        // A -> RA
                        newRA = binomial(A[i][k], mu);

                        // update
                        Snext[i][k] = S[i][k] - newL;
                        Lnext[i][k] = L[i][k] + newL - newP;
                        Pnext[i][k] = P[i][k] + newP - newA - newI;
                        Inext[i][k] = I[i][k] + newI - newRI;
                        Anext[i][k] = A[i][k] + newA - newRA;
                        RInext[i][k] = RI[i][k] + newRI;
                        RAnext[i][k] = RA[i][k] + newRA;
                    }
                }

                // update next step data and write results
                for (int i = 0; i < Npop; i++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        S[i][k]  = Snext[i][k];
                        L[i][k]  = Lnext[i][k];
                        P[i][k]  = Pnext[i][k];
                        I[i][k]  = Inext[i][k];
                        A[i][k]  = Anext[i][k];
                        RI[i][k] = RInext[i][k];
                        RA[i][k] = RAnext[i][k];
                        resFile << P[i][k] << "," << I[i][k] << "," << A[i][k] << "," << RI[i][k] << "," << RA[i][k] << ",";
                    }
                }
                resFile << "\n";
            }
        }
    }
    return 0;
}
