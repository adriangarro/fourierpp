//
//  main.cpp
//  Fourierpp
//
//  Test Discrete Fourier Transform algorithms
//  with random sine wave points as sample.
//
//  Created by Elberth Adrián Garro Sánchez on 30/3/16.
//  Copyright © 2016 Elberth Adrián Garro Sánchez. All rights reserved.
//

#include <cmath>
#include <limits>
#include <chrono>
#include <random>
#include <vector>
#include <complex>
#include <fstream>
#include <numeric>

#include "fourier.hpp"

using namespace std;
using namespace fourierpp;

using cmplx = complex<double>;
using signal_vector = vector<cmplx>;

class Fourier_analyzer {
    private:
        int menu_choice;
        signal_vector data;
        signal_vector result;
        vector<size_t> sample_sizes;
        vector<double> runtimes;
        double average_runtime;
        double variance_runtime;
        vector<string> averages_info;
        vector<string> variances_info;
        
        void validate_menu_choice()
        {
            while (cin.fail()) {
                cin.clear();
                cin.ignore(numeric_limits<int>::max(), '\n');
                cout << "Input error, please try again: ";
                cin >> menu_choice;
                cout << endl;
            }
        }
        
        void set_sample_sizes()
        {
            sample_sizes = {128, 256, 512, 1024, 2048, 4096, 8192};
        }
        
        void signal_sample(double sampling_freq, double amp1,
                           double freq1, double amp2, double freq2)
        {
            // create the points to sample in (0, 1) single period
            double time_interval = 1.0 / sampling_freq;
            vector<double> time_points;
            for (auto time_point = 0.0; time_point < 1.0; time_point += time_interval) {
                time_points.push_back(time_point);
            }
            // sum two freq Hz sine waves with peak amplitude amp (phase = 0)
            for (auto& time_point : time_points) {
                // signal(t) = sin(amp1*2*pi*freq1*t) + sin(amp2*2*pi*freq2*t)
                auto signal1 = amp1 * sin(2 * pi * freq1 * time_point);
                auto signal2 = amp2 * sin(2 * pi * freq2 * time_point);
                data.push_back(signal1 + signal2);
            }
        }
        
        int rand_int(int min, int max)
        {
            // random engine
            static random_device rd;
            return uniform_int_distribution<int>{min, max}(rd);
        }
        
        void set_data(double sample_freq)
        {
            // random signal amplitudes greater than one
            auto amp1 = rand_int(1, numeric_limits<int>::max());
            auto amp2 = rand_int(1, numeric_limits<int>::max());
            // random signal frequencies greater than one
            auto freq1 = rand_int(1, numeric_limits<int>::max());
            auto freq2 = rand_int(1, numeric_limits<int>::max());
            // create random signal on data
            signal_sample(sample_freq, amp1, freq1, amp2, freq2);
        }
        
        void verify_sdft_idft()
            // check dft with idft
        {
            cout << "\nVerification of DFT and IDFT\n";
            // create signal of 8 points on data
            set_data(8);
            cout_signal_vector("Random signal:", data);
            result = sdft(data);
            cout_signal_vector("Discrete Fourier Transform:", result);
            data = result;
            result = idft(data);
            cout_signal_vector("Inverse Discrete Fourier Transform:", result);
            data.clear();
        }
        
        void verify_fft_idft()
            // check fft with idft
        {
            cout << "\nVerification of FFT and IDFT\n";
            // create signal of 8 points on data
            set_data(8);
            cout_signal_vector("Random signal:", data);
            result = fft(data);
            cout_signal_vector("Discrete Fourier Transform:", result);
            data = result;
            result = idft(data);
            cout_signal_vector("Inverse Discrete Fourier Transform:", result);
            data.clear();
        }
        
        void save_log(string algorithm_name, string experiment_info)
        {
            ofstream log_file;
            log_file.open(algorithm_name + ".log", ios::app);
            log_file << experiment_info + "\n";
            log_file.close();
        }
        
    
        double measure_time(string algorithm_name)
            // measure dft algorithm's time by its name
        {
            chrono::high_resolution_clock::time_point time1;
            chrono::high_resolution_clock::time_point time2;
            // measure simple discrete fourier fransform algorithm
            if (algorithm_name == "sdft") {
                time1 = chrono::high_resolution_clock::now();
                result = sdft(data);
                time2 = chrono::high_resolution_clock::now();
            }
            // measure fast fourier fransform algorithm
            else if (algorithm_name == "fft") {
                time1 = chrono::high_resolution_clock::now();
                result = fft(data);
                time2 = chrono::high_resolution_clock::now();
            }
            return chrono::duration_cast<chrono::milliseconds>(time2-time1).count();
        }
        
        void set_average_runtime()
        {
            // compute unbiased average
            // E(X) = (sum from 1 to N of x) / N
            average_runtime = accumulate(runtimes.begin(), runtimes.end(), 0.0);
            average_runtime /= runtimes.size();
        }
        
        void save_average_info(size_t sample_size)
        {
            string average_info = string("Average of execution ")
                + string("time with ")
                + string("a sample of ")
                + to_string(sample_size)
                + string(" points: ")
                + to_string(int(average_runtime))
                + string(" miliseconds.");
            averages_info.push_back(average_info);
        }
        
        void set_variance_runtime()
        {
            // clear variance
            variance_runtime = 0.0;
            // compute unbiased variance
            // var(X) = [ sum from 1 to N of (x - E(X))^2 ] / N-1
            for (auto& runtime : runtimes) {
                variance_runtime += pow(runtime - average_runtime, 2);
            }
            variance_runtime /= runtimes.size() - 1;
        }
        
        void save_variance_info(size_t sample_size)
        {
            string variance_info = string("Variance of execution ")
                + string("time with ")
                + string("a sample of ")
                + to_string(sample_size)
                + string(" points: ")
                + to_string(int(variance_runtime))
                + string(" miliseconds.");
            variances_info.push_back(variance_info);
        }
        
        void do_experiments(size_t sample_size, string algorithm_name)
        {
            //-----------------------------------------------------------------
            // do one hundred experiments
            for (auto experiment_num = 1; experiment_num <= 100; ++experiment_num) {
                //-----------------------------------------------------------------
                // create signal on data
                set_data(sample_size);
                //-----------------------------------------------------------------
                // create current experiment info
                string current_experiment_info = string("Experiment number ")
                    + to_string(experiment_num)
                    + string(" with a sample of ")
                    + to_string(int(sample_size))
                    + string(" points. ");
                // record current experiment info
                save_log(algorithm_name, current_experiment_info);
                cout << current_experiment_info;
                //-----------------------------------------------------------------
                // calculate algorithm runtime
                double algorithm_runtime = measure_time(algorithm_name);
                // save runtime
                runtimes.push_back(algorithm_runtime);
                // create algorithm runtime info
                string algorithm_runtime_info = string("Execution time: ")
                    + to_string(int(algorithm_runtime))
                    + string(" miliseconds.");
                // record algorithm runtime info
                save_log(algorithm_name, algorithm_runtime_info);
                cout << algorithm_runtime_info
                     << endl;
                //-----------------------------------------------------------------
                // clear random input signal
                data.clear();
                //-----------------------------------------------------------------
            }
            // after one hundred experiments
            //-----------------------------------------------------------------
            // establish average_runtime
            set_average_runtime();
            save_average_info(sample_size);
            //-----------------------------------------------------------------
            // establish variance_runtime
            set_variance_runtime();
            save_variance_info(sample_size);
            //-----------------------------------------------------------------
            // clear runtimes
            runtimes.clear();
            //-----------------------------------------------------------------
        }
        
        void fourier_test(string algorithm_name)
        {
            //-----------------------------------------------------------------
            // establish sample sizes
            set_sample_sizes();
            //-----------------------------------------------------------------
            // create initial algorithm test info
            string initial_msg = "Algorithm's name: "
                + algorithm_name
                + "\nStarting tests...";
            // save initial algorithm test info
            save_log(algorithm_name, initial_msg);
            cout << initial_msg
                 << endl;
            //-----------------------------------------------------------------
            // separate info
            cout << endl;
            //-----------------------------------------------------------------
            // start algorithm test with point sample sizes
            for (auto sample_size : sample_sizes) {
                do_experiments(sample_size, algorithm_name);
                // separate info
                save_log(algorithm_name, "\n");
                cout << endl;
            }
            //-----------------------------------------------------------------
            // save final average results
            for (auto average : averages_info) {
                save_log(algorithm_name, average);
                cout << average
                     << endl;
            }
            //-----------------------------------------------------------------
            // separate info
            cout << endl;
            save_log(algorithm_name, "\n");
            //-----------------------------------------------------------------
            // save final variance results
            for (auto variance : variances_info) {
                save_log(algorithm_name, variance);
                cout << variance
                     << endl;
            }
            //-----------------------------------------------------------------
            // separate info
            cout << endl;
            //-----------------------------------------------------------------
            // end process
            cout << "Analysis done. The log file has been created.\n";
            averages_info.clear();
            variances_info.clear();
            //-----------------------------------------------------------------
        }
    
    public:
        Fourier_analyzer() = default;
        
        void start()
        {
            cout << endl
            << "---------------------(Menu)------------------------"
            << endl
            << "--------------Algorithm's Analysis---------------"
            << endl
            << endl
            << "            1) Verification: DFT - IDFT"
            << endl
            << "            2) Verification: FFT - IDFT"
            << endl
            << "            3) Perfomance test: Simple DTF"
            << endl
            << "            4) Perfomance test: FFT"
            << endl
            << "            0) End"
            << endl
            << "---------------------------------------------------"
            << endl
            << "Option: ";
            cin >> menu_choice;
            validate_menu_choice();
            switch (menu_choice) {
                case 1:
                    verify_sdft_idft();
                    start();
                    break;
                case 2:
                    verify_fft_idft();
                    start();
                    break;
                case 3:
                    fourier_test("sdft");
                    start();
                    break;
                case 4:
                    fourier_test("fft");
                    start();
                    break;
                case 0:
                    cout << "Thank you for your time :)"
                         << endl;
                    break;
                default:
                    cout << "Pick a valid option from menu."
                         << endl;
                    start();
            }
        }
};

Fourier_analyzer get_analyzer()
{
    return Fourier_analyzer{};
}

int main()
{
    auto analyzer = get_analyzer();
    analyzer.start();
}
