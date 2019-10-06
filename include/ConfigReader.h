//
// Created by Suman Kalyan Bera on 2019-09-27.
//

#ifndef SUBGRAPHCOUNT_CONFIGREADER_H
#define SUBGRAPHCOUNT_CONFIGREADER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;
struct config_params{
    vector<string> input_files;
    vector <long long int> triangle_count;
    int no_of_repeats = 1;
    vector <double> sparsification_prob;
    double subsample_prob;
    vector<int> seed_count;
    vector <string> algo_names;
    bool print_to_console =true;
    bool print_to_file=false;
    bool degree_bin_seed=false;
};

void ParseToken(config_params& cfp, string key, string val) {

    if (key == "input_files") {
        stringstream is_val (val);
        string token;
        while (getline(is_val,token,','))
            cfp.input_files.emplace_back(token);
    }

    else if (key == "triangle_count") {
        stringstream is_val (val);
        string token;
        while (getline(is_val,token,','))
            cfp.triangle_count.emplace_back(stod(token));
    }

    else if (key == "no_of_repeats")
        cfp.no_of_repeats = stoi(val, nullptr,10);

    else if (key == "sparsification_prob") {
        stringstream is_val (val);
        string token;
        while (getline(is_val,token,','))
            cfp.sparsification_prob.emplace_back(stod(token));
    }

    else if (key == "subsample_prob")
        cfp.subsample_prob = stod(val);

    else if ( key == "seed_count") {
        stringstream is_val (val);
        string token;
        while (getline(is_val,token,','))
            cfp.seed_count.emplace_back(stoi(token));
    }


    else if (key == "algo_names") {
        stringstream is_val (val);
        string token;
        while (getline(is_val,token,','))
            cfp.algo_names.emplace_back(token);
    }

    else if (key == "print_to_console") {
        if (val == "true")
            cfp.print_to_console = true;
        else if (val == "false")
            cfp.print_to_console = false;
        else
            std::cout << "Unknow token encountered \n";
    }
    else if (key == "print_to_file") {
        if (val == "true")
            cfp.print_to_file = true;
        else if (val == "false")
            cfp.print_to_file = false;
        else {
            std::cout << "Unknow token encountered \n" ;
        }
    }
    else if (key == "degree_bin_seed") {
        if (val == "true")
            cfp.degree_bin_seed = true;
        else if (val == "false")
            cfp.degree_bin_seed = false;
        else {
            std::cout << "Unknow token encountered \n" ;
        }
    }

    else
        std::cout << "Unknow token encountered \n" ;
}

config_params LoadConfig(const char* argv) {
    config_params cfg;
    std::ifstream ifs(argv);
    std::string line;
    if (ifs.is_open()) {
        while (std::getline(ifs,line)) {
            if (line.length() == 0 || line[0]=='#')  // Ignore empty lines and comments
                continue;
            istringstream is_line(line);
            string key,val;
            getline(is_line,key,'=');
            getline(is_line,val,'=');
            ParseToken(cfg, key, val);
        }
        ifs.close();
    }
    else {
        std::cout << "Unable to open config file \n";
        exit (-1);
    }
    return cfg;
}
#endif //SUBGRAPHCOUNT_CONFIGREADER_H
