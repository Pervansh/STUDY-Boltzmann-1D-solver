#pragma once

#include<vector>
#include<ostream>
#include<fstream>
#include<concepts>
#include<string>
#include<filesystem>

#include <BoltzmannEq1DSolvers.h>
#include <VectorOperations.h>

template<std::floating_point T>
struct Macroparameters1d;

class Macroparameters1dOutput {
private:
    std::string folderName;

    std::ofstream t_output;
    std::ofstream n_output;
    std::ofstream u_1_output;
    std::ofstream T_output;
    std::ofstream q_1_output;

public:
    Macroparameters1dOutput(const std::string& folderName) : folderName(folderName)
    {
        if (!std::filesystem::exists(folderName)) {
            if (!std::filesystem::create_directory(folderName)) {
                std::cerr << "folder \"" << folderName << "\" wasn't created!";
            }
        }
        t_output.open(folderName + "\\t_levels.txt");
        n_output.open(folderName + "\\n.txt");
        u_1_output.open(folderName + "\\u_1.txt");
        T_output.open(folderName + "\\T.txt");
        q_1_output.open(folderName + "\\q_1.txt");
    }

    template<std::floating_point T>
    void print(T t, const Macroparameters1d<T>& params) {
        t_output << t << ' ';
        t_output.flush();

        printVector(params.n, n_output);
        n_output << std::endl;

        printVector(params.u_1, u_1_output);
        u_1_output << std::endl;

        printVector(params.T, T_output);
        T_output << std::endl;

        printVector(params.q_1, q_1_output);
        q_1_output << std::endl;
    }

    void close() {
        t_output.close();
        n_output.close();
        u_1_output.close();
        T_output.close();
        q_1_output.close();
    }
};
