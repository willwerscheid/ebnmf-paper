% TODO: Explain here what this script is for, and how to use it.

% Initialize the pseudorandom number sequence.
rng(1);

% Load the swimmer data set.
load('../data/swimmer.mat');
n = 1024;
m = 256;
Y = reshape(Y,n,m);
Y = Y - 1;

% First run "vanilla" NMF without any sparsity penalty.
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
save('swimmer_nmf.mat','W','H','sp');

% k = 17; % Number of factors.
% sW = 0.95; % Average sparsity of W.

return

disp(sum(sum(W > 1e-6))); % Sparsity of W.
fprintf('Sparsity in columns of H:\n');
fprintf('%d ',sum(H > 0.001));
fprintf('\n');
clf;
sp = zeros(1,k);
for i = 1:k
  subplot(5,4,i);
  imshow(reshape(W(:,i),32,32));
  sp(i) = sp_col(W(:,i));
end
fprintf('Hoyer sparsity of columns of W: %0.3f\n',mean(sp));
