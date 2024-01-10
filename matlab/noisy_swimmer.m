% Initialize the pseudorandom number sequence.
rng(1);

% Load the noisy swimmer data set.
load('../data/noisy_swimmer.mat');

% First run "vanilla" NMF without any sparsity constraint.
k = 16;
options.sW = 0;
options.sH = 0;
options.maxiter = 100; 
options.delta = 1e-8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('noisy_swimmer_nmf.mat','W','H','sp');

% Now run NMF with a constraint on the sparsity of W. Note that since the
% unconstrained estimate of W is already quite sparse, we have to make the
% desired sparsity quite high to have any effect on the NMF estimates.
k = 17;
options.sW = 0.8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('noisy_swimmer_nmf_sW=0.8.mat','W','H','sp');

% Now run NMF with a constraint on the sparsity of W. Note that since the
% unconstrained estimate of W is already quite sparse, we have to make the
% desired sparsity quite high to have any effect on the NMF estimates.
options.sW = 0.9;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('noisy_swimmer_nmf_sW=0.9.mat','W','H','sp');

