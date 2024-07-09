%% init

% Load the noisy swimmer data set.
load('simdata2.mat');

%% vanilla NMF

% Initialize the pseudorandom number sequence.
rng(1);

% First run "vanilla" NMF without any sparsity constraint.
k = 12;
options.sW = 0;
options.sH = 0;
options.maxiter = 100; 
options.delta = 1e-8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparseness of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata2_nmf_vanilla.mat','W','H','sp');

%% Sparse NMF, sW and sH constrained

% Initialize the pseudorandom number sequence.
rng(1);

options.sW = 0.5;
options.sH = 0.7;
[W,H,e,t] = sparseNMF(Y,k,options);

sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

save('simdata2_nmf_sW=0.5.mat','W','H','sp');

options.sW = 0.6;
options.sH = 0.8;
[W,H,e,t] = sparseNMF(Y,k,options);

sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

save('simdata2_nmf_sW=0.6.mat','W','H','sp');

%% sparse NMF, greedy init

% Run vanilla NMF with the "greedy" initialization from
% flash_greedy_init_default() in R.
k = 6;
options.W = importdata('simdata_flash_greedy_init_L.txt');
options.H = importdata('simdata_flash_greedy_init_F.txt')';

options.sW = 0.4;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata_nmf_greedy_init_sW=0.4.mat','W','H','sp');

options.sW = 0.6;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata_nmf_greedy_init_sW=0.6.mat','W','H','sp');

options.sW = 0.7;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata_nmf_greedy_init_sW=0.7.mat','W','H','sp');

options.sW = 0.8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata_nmf_greedy_init_sW=0.8.mat','W','H','sp');

%% sparse NMF, both constrained

options.sW = 0.6;
options.sH = 0.8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('simdata_nmf_greedy_init_sW=0.6_sH=0.8.mat','W','H','sp');
