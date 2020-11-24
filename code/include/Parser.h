//
// Created by Staff on 22/06/2020.
//

#ifndef SIMPLE_PARSER_H
#define SIMPLE_PARSER_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "json.hpp"

using namespace std;
//using json = nlohmann::json;


/**
 * Parser Class
 */
class Parser
{
    public:
        string file;
        nlohmann::json jf;

    Parser(string f)
    {
        file = f;
        ifstream ifs(file);
        jf = nlohmann::json::parse(ifs);
    }


    /**
     * Parse commuting networks
     */
    vector<vector<double>> parse_commuting()
    {
        string commuting_file = jf["commuting_file"];
        ifstream ifs(commuting_file);
        nlohmann::json jf_commuting = nlohmann::json::parse(ifs);
        int N = jf_commuting["N"];
        vector<vector<double>> sigmas(N, vector<double>(N, 0.0));
        for (auto& el : jf_commuting["network"].items())
            sigmas[int(el.value()["from"])][int(el.value()["to"])] = double(el.value()["weight"]);
        return sigmas;
    }


    /**
     * Parse age groups contacts
     */
    vector<vector<double>> parse_contacts(int m)
    {
        string contacts_file;
        if (m==1)
            contacts_file = jf["contacts_file_home"];
        else
            contacts_file = jf["contacts_file_other"];
        ifstream ifs(contacts_file);
        nlohmann::json jf_contacts = nlohmann::json::parse(ifs);
        int K = jf_contacts["K"];
        vector<vector<double>> C(K, vector<double>(K, 0.0));
        for (auto& el : jf_contacts["rates"].items())
            C[int(el.value()["from"])][int(el.value()["to"])] = double(el.value()["rate"]);
        return C;
    }


    /**
     * Parse compartment
     */
    vector<vector<double>> parse_compartments(string c)
    {
        vector<vector<double>> compartment;
        vector<double> compart_age;
        string compartments_file = jf["populations_file"];
        ifstream ifs(compartments_file);
        nlohmann::json jf_compartments = nlohmann::json::parse(ifs);
        for (auto& el : jf_compartments["populations"].items())
        {
            compart_age.clear();
            for (auto &age_el : el.value()["age"].items())
                compart_age.push_back(age_el.value()[c]);
            compartment.push_back(compart_age);
        }
        return compartment;
    }

    /**
     * Parse r
     */
    vector<double> parse_r(int i)
    {
        vector<double> r;
        string rr = "r1";
        if (i != 1)
            rr = "r2";
        string r_file = jf["populations_file"];
        ifstream ifs(r_file);
        nlohmann::json jf_r = nlohmann::json::parse(ifs);
        for (auto& el : jf_r["populations"].items())
            r.push_back(el.value()[rr]);
        return r;
    }
};

#endif //SIMPLE_PARSER_H
