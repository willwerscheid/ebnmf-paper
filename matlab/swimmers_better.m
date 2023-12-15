rng(1);
sW = 0;
k = 16;
load('../data/swimmer.mat');
Y = reshape(Y,1024,256);
Y = Y - 1;
clear options
options.sW = sW; 
options.sH = 0;
options.maxiter = 100; 
options.delta = 1e-8;
[W0,H0] = sparseNMF(Y,k,options);
torso = 10 * (sum(W0,2) > 19);
for i = 1:k
  W0(torso > 0,i) = 0;
end
k = 17;
sW = 0;
options.sW = 0.97; 
options.sH = 0;
options.W = [W0 torso];
options.H = [H0; ones(1,256)];
[W,H,e,t] = sparseNMF(Y,k,options);
disp(sum(sum(W > 1e-6))); % Sparsity of W.
fprintf('Sparsity in columns of H:\n');
fprintf('%d ',sum(H > 0.001));
fprintf('\n');
clf;
x = zeros(1,k);
for i = 1:k
  subplot(5,4,i);
  imshow(reshape(W(:,i),32,32));
  x(i) = sp_col(W(:,i));
end
fprintf('Hoyer sparsity of columns of W: %0.3f\n',mean(x));
