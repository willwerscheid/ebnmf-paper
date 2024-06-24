% Initialize the pseudorandom number sequence.
rng(1);

% Load the noisy swimmer data set.
load('../data/noisy_swimmer.mat');

% First run "vanilla" NMF without any sparsity constraint.
k = 17;
options.sW = 0;
options.sH = 0;
options.maxiter = 100; 
options.delta = 1e-8;
[W,H,e,t] = sparseNMF(Y,k,options);

% Now run NMF with a constraint on the sparsity of W. Note that since the
% unconstrained estimate of W is already quite sparse, we have to make the
% desired sparsity quite high to have any effect on the NMF estimates.
options.sW = 0.95;
options.sH = 0.5;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('noisy_swimmer_nmf_sW=true','W','H','sp');

options.alpha = ones(1024,17);
options.beta = ones(17,256);
options.theta = 1;
options.k = 1;
[Am,Bm,sm] = gibbsnmf(Y,17,1000,options);
save('noisy_swimmer_bayesnmf','Am','Bm','sm');
