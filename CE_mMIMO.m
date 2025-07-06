clear;
close all;
clc;
format long;
%--------------------------------
% variables that will remain constant in simulations
% channel setup (Monte carlo loops)
max_avg = 50;
% Select what you want to vary with NMSE
SNR = 1;

if SNR == 1
    snr_db0 = 0: 5: 30;
    loop = length(snr_db0);
    N0 = 128;  % No. of antennas base station have
    M0 = 1;  % No. of antennas users have
    T0 = 40;  % Time slots
    K0 = 10;   % No of users
    U0 = 2;  % user specific
    C0 = 2;  % common 
    S0 = 2;  % cluster
    L0 = 2;  % No of clusters
end

% initialization of error metrics
nmse = zeros(loop,1);
%----------------------------------------
%========================================
% select the algorithms you want to run
channel_est = 1;
%========================================
%----------------------------------------
% Progress bar
wt = waitbar(0,'Initializing...');
curr_iter = 0;
for loop_iter = 1:loop
    nmse_sum = 0;
    for mc_iter = 1:max_avg
        if SNR == 1
            snr = 10.^(snr_db0(loop_iter)/10);
            N = N0;
            M = M0;
            T = T0;
            K = K0;
            U = U0;
            C = C0;
            S = S0;
            L = L0;
        end
        % X_true
        CL = create_clusters(K, L);
        A_R = dftmtx(M)/sqrt(M);
        A_T = dftmtx(N)/sqrt(N);
        [X_true, H] = Channel_Generation(N, M, K, U, C, S, CL, A_R, A_T);
        sigma2_k = 1/snr;
        % Phi = sqrt(1/N)*(unidrnd(2,[T,N])-1.5)*2;
        Phi = randn(T,N);
        % Phi = (1/sqrt(N))*exp(1i*2*pi*rand(T,N));
        Y = zeros(T, M, K);
        for k = 1 : K
            Xk = X_true(:,:,k);
            Yk = signal_gen_mMIMO(Phi, Xk, M, T, sigma2_k, A_R);
            Y(:, :, k) = Yk;
        end

        if channel_est == 1
            X_est = updates(Phi, Y, sigma2_k, CL,X_true);
            H_est = zeros(M, N,K);
            for k = 1 : K
                H_est(:,:,k) = A_R*(X_est(:,:,k)')*A_T';
            end
            for k = 1 : K
                nmse_sum = nmse_sum + (norm(X_est(:, :, k) - X_true(:, :, k),'fro')^2 / norm(X_true(:, :, k),'fro')^2)/K;
            end
        end
        %% ----------------------------------------
        % Update waitbar
        curr_iter = 1 + curr_iter;
        waitbar(curr_iter/(max_avg*loop),wt,sprintf('%0.1f%% done',curr_iter/(max_avg*loop)*100))
    end
    nmse(loop_iter) = nmse_sum / max_avg;

    % Plot result
    if SNR == 1
        semilogy(snr_db0, nmse, '-o', 'LineWidth', 2, 'DisplayName', 'ChannelEstimation')
        xlabel('SNR (dB)');
    end
    ylabel('NMSE');
    legend show
    grid on
end
delete(wt)