/* 
    Implements Discrete Fourier Transform Algorithms
    
    Built with C++14
    
    (c) E. Adrian Garro S. Costa Rica Institute of Technology.
 */

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>

namespace fourierpp {
    static constexpr double pi {3.1415926535897932384626433832795};

    template <class signal_vector>
    void cout_signal_vector(std::string msg, signal_vector& data)
    {
        std::cout << msg
                  << std::endl; 
        for (auto point : data) {
            std::cout << point
                      << std::endl;
        }
    }
    
    void adjust_zero_precision(std::complex<double>& point)  
        // close to zero? you're zero
    {
      auto real_part = point.real();
      auto imag_part = point.imag();
      if (std::abs(real_part) < 1e-5) {
        real_part = 0.0;
        point = std::complex<double>(real_part, imag_part);
      }
      if (std::abs(imag_part) < 1e-5) {
        imag_part = 0.0;
        point = std::complex<double>(real_part, imag_part);
      }
    }
    
    template <class signal_vector>
    void adjust_result(signal_vector& result)
        // points of signal close to zero become zero 
    {
        for (auto& point : result) {
            adjust_zero_precision(point);
        }
    }
    
    template <class signal_vector>
    signal_vector sdft(signal_vector& data) 
        // Simple DFT (not in-place)
    {
        auto data_size = data.size();
        signal_vector result(data_size);
        // for every frequency...
        for (auto freq = 0; freq < data_size; ++freq) {
            // for every point in time...
            for (auto t = 0; t < data_size; ++t) {  
                auto twiddle_factor = std::polar(1.0, -2 * pi * freq * t / data_size);
                // datapoint * e^(-i*2*pi*f)
                auto result_point = data[t] * twiddle_factor;
                result[freq] += result_point;
            }
        }
        adjust_result(result);
        return result;
    }
    
    template <class signal_vector, typename size>
    void partition_by_index(signal_vector& data, size data_size, 
                            signal_vector& odd, signal_vector& even) 
    {
        for (auto index = 0; index < data_size; ++index) {
            if ((index % 2) != 0) odd.push_back(data.at(index));
            else even.push_back(data.at(index));
        }
    }
    
    template <class signal_vector>
    void aux_fft(signal_vector& data)  
        // Cooley–Tukey FFT (in-place)
    {
        auto data_size = data.size();
        if (data_size == 1) return;
        // divide
        signal_vector odd, even;
        partition_by_index(data, data_size, odd, even);
        // conquer
        aux_fft(even);
        aux_fft(odd);
        // combine
        // for every frequency...
        for (auto freq = 0; freq < data_size / 2; ++freq) {
            auto twiddle_factor = std::polar(1.0, -2.0 * pi * freq / data_size);
            data[freq] = even[freq] + twiddle_factor * odd[freq];
            data[freq + data_size / 2] = even[freq] - twiddle_factor * odd[freq];
        }
    }
    
    template <class signal_vector>
    signal_vector fft(signal_vector& data)  
        // Cooley–Tukey FFT (not in-place)
    {
        signal_vector result(data);
        aux_fft(result);
        return result;
    }
    
    template <class signal_vector>
    signal_vector idft(signal_vector& data)  
        // Inverse DFT (not in-place)
    {
        auto data_size = data.size();
        signal_vector result(data_size);
        // for every point in time...
        for (auto t = 0; t < data_size; ++t) {
            // for every frequency...
            for (auto freq = 0; freq < data_size; ++freq) {
                auto twiddle_factor = std::polar(1.0, 2.0 * pi * freq * t / data_size);
                // datapoint * e^(i*2*pi*f*t/N)
                auto result_point = data[freq] * twiddle_factor;
                result[t] += result_point;
            }
            result[t] /= data_size;
        }
        adjust_result(result);
        return result;
    }
}
