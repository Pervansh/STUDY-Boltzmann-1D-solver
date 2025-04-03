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

template<std::floating_point T>
struct Destribution1dState;

class Macroparameters1dOutput {
private:
    std::string folderName;

    std::ofstream t_output;
    std::ofstream n_output;
    std::ofstream u_1_output;
    std::ofstream T_output;
    std::ofstream q_1_output;
    std::ofstream x_1_output;

public:
    Macroparameters1dOutput(const std::string& folderName) : folderName(folderName)
    {
        if (!std::filesystem::exists(folderName)) {
            if (!std::filesystem::create_directory(folderName)) {
                std::cerr << "folder \"" << folderName << "\" wasn't created!";
            }
        }

        t_output.open(folderName + "\\t_levels.txt");
        x_1_output.open(folderName + "\\x_1_levels.txt");

        n_output.open(folderName + "\\n.txt");
        u_1_output.open(folderName + "\\u_1.txt");
        T_output.open(folderName + "\\T.txt");
        q_1_output.open(folderName + "\\q_1.txt");
    }

    template<std::floating_point T>
    void print(T t, const std::vector<T>& x_v, const Macroparameters1d<T>& params) {
        t_output << t << ' ';
        t_output.flush();

        printVector(x_v, x_1_output);
        x_1_output << std::endl;

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
        x_1_output.close();
        n_output.close();
        u_1_output.close();
        T_output.close();
        q_1_output.close();
    }
    
    std::string getFolderName() { return folderName; };
};

class Full1dStateOutput {
private:
    Macroparameters1dOutput mpOutput;

    std::ofstream h_output;
    std::ofstream g_output;
    std::ofstream xi_1_output;

public:
    Full1dStateOutput(const std::string& folderName) : mpOutput(folderName)
    {
        xi_1_output.open(folderName + "\\xi_1_levels.txt");

        h_output.open(folderName + "\\h.txt");
        g_output.open(folderName + "\\g.txt");
    }

    /**
    * /param params must be consistent with state
    */
    template<std::floating_point T>
    void print(
        const T t, 
        const std::vector<T>& x_v,
        const std::vector<T>& xi_v,
        const Macroparameters1d<T>& params,
        const Destribution1dState<T>& state
    ) {
        mpOutput.print<T>(t, x_v, params);

        printVector(xi_v, xi_1_output);
        xi_1_output << std::endl;

        printVector(state.h, h_output);
        h_output << std::endl << std::endl;

        printVector(state.g, g_output);
        g_output << std::endl << std::endl;
    }

    void close() {
        mpOutput.close();

        xi_1_output.close();
        h_output.close();
        g_output.close();
    }

    std::string getFolderName() { return mpOutput.getFolderName(); };
};
