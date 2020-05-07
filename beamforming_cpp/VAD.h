#include <iostream>
#include <math.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/MatrixFunctions>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

//--------------------------- Hamming windows ------------------------------//

MatrixXd hamming(int const n)
{
    MatrixXd t(n, 1);

    if (n == 1)
        t(0, 0) = 0.08;
    else
        for (int i = 0; i < n; i++)
            t(i, 0) = (0.54 - 0.46 * cos(2.0 * M_PI * i / (n - 1)));

    return t;
}

//--------------------------- Besseli: modified Bessel function ------------------------------//

MatrixXd besseli(double x, MatrixXd v)
{
    int rows = v.rows(), cols = v.cols();
    MatrixXd result(rows, cols);
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            result(r, c) = cyl_bessel_i(x, v(r, c));
        }
    }
    return result;
}

//--------------------------- VAD ------------------------------//

MatrixXi VAD(MatrixXd &audio, int fs, double threshold, double win_dur, double hop_dur, int num_noise, int argin)
{
    int DEBUG = 1;

    // Markov parameters for hangover scheme
    double a01 = 0.5;
    double a10 = 0.1;
    double a00 = 1 - a01;
    double a11 = 1 - a10;

    // Coefficient for decision-directed SNR estimation
    double alpha = 0.9900;

    int a = audio.cols();
    int audio_lenght = audio.rows();
    int win_size = round(fs * win_dur);
    int hop_size = round(fs * hop_dur);
    int N_FFT = win_size;
    int num_frames = floor((audio_lenght - win_size) / hop_size);
    FFT<double> fft;

    // Const matrix
    MatrixXd one(N_FFT, 1), zero(N_FFT, 1), alpha_matrix(N_FFT, 1);
    one.setOnes();
    zero.setZero();
    alpha_matrix.fill(alpha);

    MatrixXd hamming_win = hamming(win_size);
    MatrixXi results = MatrixXi::Zero(num_frames, 1);
    MatrixXd noise_var_temp = MatrixXd::Zero(N_FFT, 1).adjoint();

    for (int i = 0; i < num_noise; i++)
    {
        MatrixXd audio_frame = hamming_win.cwiseProduct(audio.block(i * hop_size, 0, win_size, 1)); // size (win_size, 1)

        VectorXcd audio_frame_fft(win_size);
        fft.fwd(audio_frame_fft, audio_frame.col(0), N_FFT);
        audio_frame_fft = audio_frame_fft.adjoint();
        VectorXcd audio_frame_fft_ctranspose = audio_frame_fft.adjoint();

        noise_var_temp.row(0) = noise_var_temp.row(0) + audio_frame_fft.cwiseProduct(audio_frame_fft_ctranspose).real().transpose();

        results(i) = 0;
    }

    MatrixXd noise_var_orig = noise_var_temp.adjoint() / num_noise;
    MatrixXd noise_var_old = noise_var_orig;

    double G_old = 1;
    MatrixXd A_MMSE = MatrixXd::Zero(N_FFT, num_frames);
    MatrixXd G_MMSE = MatrixXd::Zero(N_FFT, num_frames);

    MatrixXd cumulative_Lambda = MatrixXd::Zero(num_frames, 1);

    for (int i = 0; i < num_frames; i++)
    {

        MatrixXd audio_frame = hamming_win.cwiseProduct(audio.block(i * hop_size, 0, win_size, 1));

        VectorXcd audio_frame_fft(win_size);
        fft.fwd(audio_frame_fft, audio_frame.col(0), N_FFT);
        audio_frame_fft = audio_frame_fft.adjoint();
        VectorXcd audio_frame_fft_ctranspose = audio_frame_fft.adjoint();

        VectorXd frame_var = audio_frame_fft.cwiseProduct(audio_frame_fft_ctranspose).real().transpose();
        MatrixXd noise_var = noise_var_orig;

        VectorXd Y_mag(N_FFT);
        MatrixXd xi(N_FFT, 1);

        if (argin)
        { // for noise estimation, set argin = 1, otherwise 0.
            int nn = 1;
            MatrixXd noise_var_prev = noise_var_orig;
            while (nn <= 10)
            {
                MatrixXd gamma = frame_var.cwiseQuotient(noise_var);
                Y_mag = audio_frame_fft.cwiseAbs();

                if (i == 0)
                {
                    // MATLAB: xi = alpha + (1-alpha) * max(gamma-1, 0);
                    xi = alpha_matrix + (one - alpha_matrix).cwiseProduct((gamma - one).cwiseMax(zero));
                }
                else
                {
                    // MATLAB: xi = alpha * ( (A_MMSE(:,i-1).^2) ./ noise_var_old ) + (1-alpha) * max(gamma-1, 0);
                    MatrixXd A_MMSE_square = A_MMSE.col(i - 1).array().square();
                    xi = alpha_matrix.cwiseProduct(A_MMSE_square.cwiseQuotient(noise_var_old)) + (one - alpha_matrix).cwiseProduct((gamma - one).cwiseMax(zero));
                }

                // MATLAB: v = (xi .* gamma) ./ ( 1 + xi);
                MatrixXd v = (xi.cwiseProduct(gamma)).cwiseQuotient(one + xi);

                // MATLAB: G_MMSE(:,i) = (sqrt(pi)/2) * (sqrt(v)./gamma) .* exp(v/-2) .*((1+v) .* besseli(0,v/2) + v .* besseli(1,v/2));
                MatrixXd exp_v = (v / -2).array().exp();
                MatrixXd besseli_part = (one + v).cwiseProduct(besseli(0, v / 2)) + v.cwiseProduct(besseli(1, v / 2));
                G_MMSE.col(i) = (sqrt(M_PI) / 2) * (v.cwiseSqrt().cwiseQuotient(gamma)).cwiseProduct(exp_v).cwiseProduct(besseli_part);

                // MATLAB: G_MMSE(find(isnan(G_MMSE(:,n))),n) = 1;
                //         G_MMSE(find(isinf(G_MMSE(:,n))),n) = 1;
                for (int k = 0; k < N_FFT; k++)
                {
                    if (isnan(G_MMSE(k, i)) || isinf(G_MMSE(k, i)))
                        G_MMSE(k, i) = 1;
                }

                // MATLAB: A_MMSE(:, n) = G_MMSE(:, n) .* Y_mag';
                A_MMSE.col(i) = G_MMSE.col(i).cwiseProduct(Y_mag.adjoint().transpose());

                MatrixXd Lambda = MatrixXd::Zero(N_FFT, 1);
                for (int k = 0; k < N_FFT; k++)
                {
                    Lambda(k, 0) = 1 / (1 + xi(k, 0)) + exp(gamma(k, 0) * xi(k, 0) / (1 + xi(k, 0)));
                }

                double Lambda_mean = Lambda.sum() / N_FFT;
                double weight = Lambda_mean / (1 + Lambda_mean);

                if (isnan(weight))
                {
                    weight = 1;
                }

                noise_var = weight * noise_var_orig + (1 - weight) * frame_var;

                double diff = abs((noise_var - noise_var_prev).sum());
                if (diff < 0.000001)
                {
                    nn = 10;
                }

                nn = nn + 1;
                noise_var_prev = noise_var;

            } //while (nn < 10)

        } // if (argin)

        MatrixXd gamma = frame_var.cwiseQuotient(noise_var);

        Y_mag = audio_frame_fft.cwiseAbs();

        if (i == 0)
        {
            // MATLAB: xi = alpha + (1-alpha) * max(gamma-1, 0);
            xi = alpha_matrix + (one - alpha_matrix).cwiseProduct((gamma - one).cwiseMax(zero));
        }
        else
        {
            // MATLAB: xi = alpha * ( (A_MMSE(:,i-1).^2) ./ noise_var_old ) + (1-alpha) * max(gamma-1, 0);
            MatrixXd A_MMSE_square = A_MMSE.col(i - 1).array().square();
            xi = alpha_matrix.cwiseProduct(A_MMSE_square.cwiseQuotient(noise_var_old)) + (one - alpha_matrix).cwiseProduct((gamma - one).cwiseMax(zero));
        }

        // MATLAB: v = (xi .* gamma) ./ ( 1 + xi);
        MatrixXd v = (xi.cwiseProduct(gamma)).cwiseQuotient(one + xi);

        // MATLAB: G_MMSE(:,i) = (sqrt(pi)/2) * (sqrt(v)./gamma) .* exp(v/-2) .*((1+v) .* besseli(0,v/2) + v .* besseli(1,v/2));
        MatrixXd exp_v = (v / -2).array().exp();
        MatrixXd besseli_part = (one + v).cwiseProduct(besseli(0, v / 2)) + v.cwiseProduct(besseli(1, v / 2));
        G_MMSE.col(i) = (sqrt(M_PI) / 2) * (v.cwiseSqrt().cwiseQuotient(gamma)).cwiseProduct(exp_v).cwiseProduct(besseli_part);

        // MATLAB: G_MMSE(find(isnan(G_MMSE(:,n))),n) = 1;
        //         G_MMSE(find(isinf(G_MMSE(:,n))),n) = 1;
        for (int k = 0; k < N_FFT; k++)
        {
            if (isnan(G_MMSE(k, i)) || isinf(G_MMSE(k, i)))
                G_MMSE(k, i) = 1;
        }

        // MATLAB: A_MMSE(:, n) = G_MMSE(:, n) .* Y_mag';
        A_MMSE.col(i) = G_MMSE.col(i).cwiseProduct(Y_mag.adjoint().transpose());

        MatrixXd Lambda = MatrixXd::Zero(N_FFT, 1);
        for (int k = 0; k < N_FFT; k++)
        {
            Lambda(k, 0) = log(1 / (1 + xi(k, 0))) + gamma(k, 0) * xi(k, 0) / (1 + xi(k, 0));
        }

        double Lambda_mean = Lambda.sum() / N_FFT;

        double G = (a01 + a11 * G_old) / (a00 + a10 * G_old) * Lambda_mean;

        MatrixXd log_gamma = gamma.array().log();
        double L_ML = (1 / N_FFT) * (gamma - log_gamma - one).sum();
        if (DEBUG)
        {
            if ((G < threshold) || (i < num_noise))
                results(i, 0) = 0;
            else
                results(i, 0) = 1;
        }

        cumulative_Lambda(i) = Lambda_mean;
        G_old = G;
        noise_var_old = noise_var;
    }

    return results;
}
