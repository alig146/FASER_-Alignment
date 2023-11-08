#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>

int Res_to_txt_5_var()
{
    std::string sensor1;
    std::string sensor2;
    double align_param1;
    double align_param2;
    std::string line;
    std::string word;
    std::string module;
    std::vector<std::string> params(6);
    std::string chars = "[]:,\"";
    std::ifstream in_file_res;
    std::ifstream in_file_input;
    std::ifstream temp_input_in;
    std::ifstream temp_res_in;
    std::ofstream temp_res_out;
    std::ofstream temp_input_out;
    std::ofstream out_file;


    in_file_res.open("/home/agarabag/millepede/target/new_alg_v2_iter1.res");
    // in_file_input.open("/home/agarabag/millepede/ift_mc_align/new_alg_iter3.txt");
    in_file_input.open("/home/agarabag/condor_jobs/ke_spoint/inputforalign.txt");
    temp_res_out.open("temporary_res.txt");
    temp_input_out.open("temporary_input.txt");
    out_file.open("new_alg_iter1_v2.txt");
    if (!in_file_res){
        std::cerr << "Could not open input millepede file" << std::endl;
        return 1;
    }
    if (!in_file_input){
        std::cerr << "Could not open input for align" << std::endl;
        return 1;
    }
    if (!out_file){
        std::cerr << "Could not open output file" << std::endl;
        return 1;
    }
    
    bool beg{true};
    
    while (std::getline(in_file_res,line)) {
        if (beg){
            beg = false;
            continue;
        } else {
            temp_res_out << line.substr(0,29) << std::endl;
        }
    }
        
    while (in_file_input >> module >> params.at(0) >> params.at(1) >> params.at(2) >> params.at(3) >> params.at(4) >> params.at(5)){
        module.erase(remove_if(module.begin(), module.end(),
            [&chars](const char &c) {return chars.find(c) != std::string::npos;}), module.end());
        int l {1};
        for (std::string param: params){
            param.erase(remove_if(param.begin(), param.end(),
                [&chars](const char &c) {return chars.find(c) != std::string::npos;}), param.end());
            // if we are at the third param skip 
            if (l == 3){
                l++;
                continue;
            }
            temp_input_out << " " << module << l << " " << param << std::endl;
            l++;
        }
    }
    
    int n {0};
    temp_res_out.close();
    temp_res_in.open("temporary_res.txt");
    temp_input_in.open("temporary_input.txt");
    while (temp_res_in >> sensor1 >> align_param1 && temp_input_in >> sensor2 >> align_param2){
        switch(sensor1[0]){
            case '1':
                sensor1[0] = '0';
                break;
            case '2':
                sensor1[0] = '1';
                break;
            case '3':
                sensor1[0] = '2';
                break;
            case '4':
                sensor1[0] = '3';
                break;
            default: break;
        }
        
        if (sensor1.length() == 3){
            if (n == 0)
                out_file << " \"" << sensor1[0] << sensor1[1] << "\": [";
            if (n < 4 && n != 1){
                 out_file << align_param1 + align_param2 << ", ";
                 n++;
            }
            else if (n == 1){
                out_file << align_param1 + align_param2 << ", 0, ";
                n++;
            } 
            else if (n == 4){
                n = 0;                
                out_file << align_param1 + align_param2 << "],";
            }
        }
        if (sensor1.length() == 5){
            std::cout << "NOOOOOOO" << std::endl;
            if (n == 0){
                out_file << " \"" << sensor1[0] << sensor1[1] << sensor1[2] << "\": [";
                out_file << align_param1 + align_param2 << ", ";
                n++;
            }
            if (n == 2 || n == 3){
                 out_file << align_param1 + align_param2 << ", ";
                 n++;
            }
            else if (n == 1){
                out_file << align_param1 + align_param2 << ", 0, ";
                n++;
            } 
            else if (n == 4){
                n = 0;
                out_file << align_param1 + align_param2 << "],";
            }
        }
    }
    temp_res_in.close();
    in_file_res.close();
    out_file.close();
    
    return 0;
}
