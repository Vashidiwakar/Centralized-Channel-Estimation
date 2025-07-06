function mean = updates(Phi, Y, sigma2_k, CL, X)
% Initialization
L = length(CL);
[~, N] = size(Phi);
[~, M, K] = size(Y);
mean  = ones(N, M, K);

alpha_k = ones(N,K)*100;
alpha_l = ones(N,L)*100;
alpha_c = ones(N,1)*100;

eta_k = zeros(N,1);
eta_l = zeros(N,1);
eta_c = zeros(N,1);

z_k = (1/3)*ones(N,1);
z_l = (1/3)*ones(N,1);
z_c = (1/3)*ones(N,1);
SIGMA = zeros(N, N, K);
mean_new = zeros(N, M, K);

max_itr = 150;
condition = 1;
i = 1;
nmse_err = zeros(max_itr,1);
stop_thres = 1e-3;
while i <= max_itr && condition == 1
    % E step : Computing mean and covariance ========
    for l = 1 : L
        for k  = CL{l}
            Yk = Y(:,:,k);
            lambda_k_diag = (z_k .* alpha_k(:,k)) + (z_l .* alpha_l(:,l)) + (z_c .* alpha_c);
            lambda_k = diag(lambda_k_diag);
            % lambda_k_inv = inv(lambda_k);
            % SIGMA_k = (lambda_k + (1/sigma2_k)*(Phi')*Phi)^(-1);
            % SIGMA_k = lambda_k_inv - lambda_k_inv * Phi' * ( (1/sigma2_k) * eye(T) + Phi * lambda_k_inv * Phi') * Phi * lambda_k_inv;
            % inv(A)*b = A\b and b*inv(A) = b/A
            % SIGMA_k = lambda_k_inv - (lambda_k \ Phi') * (((1/sigma2_k) * eye(T) + Phi *(lambda_k \ Phi'))\(Phi / lambda_k));
            A = lambda_k + (1/sigma2_k)*(Phi')*Phi;
            SIGMA_k = A \ eye(N);  % More stable than inv(A)

            % for m = 1 : M
            %     mean_new(:,m,k) = (1/sigma2_k) * SIGMA_k * (Phi') * Y(:, m, k); 
            % end
            mean_new(:,:,k) = (1/sigma2_k)*SIGMA_k*Phi'*Yk;
            SIGMA(:, :, k) = SIGMA_k;
        end
    end

    % computing eta_k, eta_l, eta_c
    
    for n = 1:N
        log_alpha_c = log(max(alpha_c(n), eps));
        eta_k(n) = 0;
        eta_l(n) = 0;
        eta_c(n) = K * M * log_alpha_c;

        for l = 1:L
            cl_size = numel(CL{l});
            for k = CL{l}  % Cluster of user k
                log_alpha_k = log(max(alpha_k(n,k), eps));
                eta_k(n) = eta_k(n) + M * log_alpha_k;
            end
            log_alpha_l = log(max(alpha_l(n,l), eps));
            eta_l(n) = eta_l(n) + M * cl_size * log_alpha_l;
        end
        % computing <z_k>, <z_l> and <z_c>
        % z_k(n,1) = 1/(1 + exp(eta_c(n,1) - eta_k(n,1)) + exp(eta_l(n,1) - eta_k(n,1)));
        % z_l(n,1) = 1/(1 + exp(eta_c(n,1) - eta_l(n,1)) + exp(eta_k(n,1) - eta_l(n,1)));
        % z_c(n,1) = 1/(1 + exp(eta_l(n,1) - eta_c(n,1)) + exp(eta_k(n,1) - eta_c(n,1)));
        den = exp(eta_k(n)) + exp(eta_l(n)) + exp(eta_c(n));
        z_k(n) = exp(eta_k(n)) / den;
        z_l(n) = exp(eta_l(n)) / den;
        z_c(n) = exp(eta_c(n)) / den;

    end

    % M step : update alpha_k, alpha_l and alpha_c ==========
    % =====Update alpha_k =====
    % for k = 1:K
    %     for n = 1:N
    %         squared_mean = 0;
    %         for m = 1:M
    %         % |mean(n,m,k)|^2 for all m
    %         squared_mean = squared_mean + abs(mean_new(n, m, k)).^2 + SIGMA(n, n, k);  % 1 x M
    %         end
    %         % Expected value = mean squared + variance
    %         alpha_k(n, k) = M / squared_mean;  % add variance
    %     end
    % end

    % =====Update alpha_k =====
    for k = 1:K
        for n = 1:N
            squared_mean = sum(abs(mean_new(n, :, k)).^2);  % mean squared term
            variance_term = M * SIGMA(n, n, k);             % variance added once for all M
            denom = squared_mean + variance_term;
            alpha_k(n, k) = M / max(denom, eps);  
        end
    end

    % =====Update alpha_l =====
    for n = 1 : N
        for l = 1 : L
            cl_size = numel(CL{l}); % gives number of element in each cluster
            alpha_k_sum = 0;
            for k = CL{l}
                alpha_k_sum = alpha_k_sum + 1/alpha_k(n,k);
            end
            alpha_l(n,l) = cl_size / max(alpha_k_sum,eps);
        end
    end
    % =====Update alpha_c =====
    for n = 1 : N
        alpha_k_sum2 = 0;
        for l = 1 : L
            for k = CL{l}
                alpha_k_sum2 = alpha_k_sum2 + 1/alpha_k(n,k);
            end
        end
        alpha_c(n,1) = K / max(alpha_k_sum2,eps);
    end
    % condition = (norm(mean_new(:) - mean(:), 2) > 1e-4);
    flag = 0;
    for k = 1 : K
        if(norm(mean_new(:,:,k) - mean(:,:,k), 2) > stop_thres*norm(mean(:,:,k), 2))
            flag = 1;
            break
        end
    end
    if(flag == 0)
        break
    end
    nmse_sum = 0;
    for k = 1 : K
        nmse_sum = nmse_sum + (norm(mean_new(:, :, k) - X(:, :, k),'fro')^2 / norm(X(:, :, k),'fro')^2)/K;
    end
    nmse_err(i) = nmse_sum;
    mean = mean_new;
    i = i + 1;
end
end