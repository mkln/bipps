#ifndef STREAMING_WAIC_H
#define STREAMING_WAIC_H

#include <armadillo>
#include <cmath>
#include <vector>
#include <limits>

class StreamingLogMeanExp {
private:
    double max_val;       // Running maximum
    double scaled_sum;    // Scaled sum of exponentials
		int count;					 // Number of samples

public:
    StreamingLogMeanExp() : max_val(-std::numeric_limits<double>::infinity()), scaled_sum(0.0), count(0) {}

    void update(double x) {
				count++;
        if (x > max_val) {
            scaled_sum = scaled_sum * std::exp(max_val - x) + 1.0;
            max_val = x;
        } else {
            scaled_sum += std::exp(x - max_val);
        }
    }

    double get_logmeanexp() const {
        return (max_val + std::log(scaled_sum / count));
    }
};

class StreamingVariance {
private:
    double mean;
    double M2;  // Sum of squared differences
    int count;

public:
    StreamingVariance() : mean(0.0), M2(0.0), count(0) {}

    void update(double x) {
        count++;
        double delta = x - mean;
        mean += delta / count;
        M2 += delta * (x - mean);
    }

    double get_variance() const {
        return count > 1 ? M2 / (count - 1) : 0.0;
    }
};

class StreamingWAIC {
private:
    int N;  // Number of observations
    std::vector<StreamingLogMeanExp> lppd_streams;
    std::vector<StreamingVariance> variance_streams;

public:
    StreamingWAIC(int num_obs) : N(num_obs) {
        lppd_streams.resize(N);
        variance_streams.resize(N);
    }

    // Update with an arma::vec of log-likelihoods (one value per observation)
    void update(const arma::vec& log_likelihood) {
        if (log_likelihood.n_elem != static_cast<arma::uword>(N)) {
            throw std::invalid_argument("Size of log_likelihood vector must match the number of observations (N).");
        }

        for (int n = 0; n < N; ++n) {
            lppd_streams[n].update(log_likelihood[n]);
            variance_streams[n].update(log_likelihood[n]);
        }
    }

    // Calculate the WAIC after all updates
    double calculate_waic() const {
        double lppd_total = 0.0;
        double penalty_total = 0.0;

        for (int n = 0; n < N; ++n) {
            lppd_total += lppd_streams[n].get_logmeanexp();
            penalty_total += variance_streams[n].get_variance();
        }

        return -2 * (lppd_total - penalty_total);
    }
};

#endif // STREAMING_WAIC_H
