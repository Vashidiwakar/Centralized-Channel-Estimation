% For unequal number of users in each cluster
function C = create_clusters(K, L)
%   C - cell array of clusters, C{1}, C{2}, ..., C{L}

% Randomly shuffle user indices
avail_idx = randperm(K);

% Compute base cluster size and remainder
base_size = floor(K / L);
extra = mod(K, L);  % extra users to distribute

C = cell(1, L);
idx = 1;

for i = 1:L
    cluster_size = base_size + (i <= extra);  % add one if extra remains
    C{i} = avail_idx(idx : idx + cluster_size - 1);
    idx = idx + cluster_size;
end
end

%=======================================================
% function C = create_clusters(K,L)
%      n = K / L; % no of users per cluster
%      C = zeros(L,n);
%      avail_idx = 1 : K;
%      for i = 1 : L
%          c_idx = avail_idx(randperm(length(avail_idx),n));
%          C(i, :) = c_idx;
%          avail_idx = setdiff(avail_idx, c_idx);
%      end
% end

