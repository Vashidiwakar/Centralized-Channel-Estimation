function [X_true,Hk] = Channel_Generation(N, M, K, U0, C0, S0, C, A_R, A_T)

    X_true = zeros(N, M, K);  % 2D array for each user, K represents the number of users
    Hk = zeros(M, N, K);
    % Common sparsity
    all_idx = 1:N;
    C_idx = all_idx(randperm(length(all_idx), C0));
    for k = 1 : K
        X_true(C_idx,:,k) = sqrt(1/2)*(randn(C0, M) + 1i * randn(C0, M));
    end
    
    % Cluster sparsity
    L = length(C);
    cluster_idx = setdiff(all_idx, C_idx);
    for i = 1 : L
        S_idx = cluster_idx(randperm(length(cluster_idx), S0));  % S_idx is numeric
        cluster_idx = setdiff(cluster_idx, S_idx);
    
        users_in_cluster = C{i};
        for k = users_in_cluster
            X_true(S_idx, :, k) = sqrt(1/2)*(randn(S0, M) + 1i * randn(S0, M));
        end
    end
    
    % User-specific sparsity
    n_idx = cluster_idx;
    
    for k = 1:K
        U_idx = n_idx(randperm(length(n_idx), U0));
        X_true(U_idx,:,k) = sqrt(1/2)*(randn(U0, M) + 1i*randn(U0, M));
        n_idx = setdiff(n_idx, U_idx);
    end
    for k = 1 : K
        Hk(:,:,k) = A_R * X_true(:,:,k)' * A_T';
    end
end