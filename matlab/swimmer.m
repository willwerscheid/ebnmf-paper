% Initialize the pseudorandom number sequence.
rng(1);

% Load the swimmer data set.
load('../data/swimmer.mat');
n = 1024;
m = 256;
Y = reshape(Y,n,m);
Y = Y - 1;

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
save('swimmer_nmf.mat','W','H','sp');

% Now run NMF with a constraint on the sparsity of W. Note that since the
% unconstrained estimate of W is already very sparse, we have to make the
% desired sparsity quite high to have any effect on the NMF estimates.
k = 17;
options.sW = 0.95;
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('swimmer_nmf_sW=0.95.mat','W','H','sp');

% Run vanilla NMF with the "greedy" initialization from
% flash_greedy_init_default() in R.
k = 17;
options.sW = 0;
options.W = importdata('flash_greedy_init_L.txt');
options.H = importdata('flash_greedy_init_F.txt')';
[W,H,e,t] = sparseNMF(Y,k,options);

% Compute the sparsity of each column of W.
sp = zeros(1,k);
for i = 1:k
  sp(i) = sp_col(W(:,i));
end

% Save the results.
save('swimmer_nmf_greedy_init.mat','W','H','sp');

%
% disp(sum(sum(W > 1e-6)));
% fprintf('Sparsity in columns of H:\n');
% fprintf('%d ',sum(H > 0.001));
% fprintf('\n');
% clf;
% sp = zeros(1,k);
% for i = 1:k
%   subplot(5,4,i);
%   imshow(reshape(W(:,i),32,32));
%   sp(i) = sp_col(W(:,i));
% end
% fprintf('Hoyer sparsity of columns of W: %0.3f\n',mean(sp));
%
