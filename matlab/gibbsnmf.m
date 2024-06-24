function [Am,Bm,sigmam] = gibbsnmf(X,N,M,varargin)
% GIBBSNMF Non-negative matrix factorization Gibbs sampler
%
% Usage
%   [Am,Bm,sm] = gibbsnmf(X,N,M,[options])
%
% Input
%   X           Data matrix (I x J)
%   N           Number of components 
%   M           Number of Gibbs samples to compute
%   options
%     alpha     Prior for A (I x N)
%     beta      Prior for B (N x J)
%     theta     Prior for sigma
%     k         Prior for sigma
%
%     A         Initial value for A (I x N)
%     B         Initial value for B (N x J)
%     sigma     Initial value for noise variance (sigma^2)
%
%     chains    Number of chains to run (default 1)
%     skip      Number of initial samples to skip (default 100)
%     stride    Return every n'th sample (default 1)
%
%     nA        Do not sample from these columns of A (N x 1 logical)
%     nB        Do not sample from these rows of B (N x 1 logical)
%     ns        Do not sample from sigma (logical)
%
% Output
%   Am          Samples of A (I x N x M)
%   Bm          Samples of B (N x J x M)
%   sm          Samples of sigma (M x 1)
%
% Author
%   Mikkel N. Schmidt, 
%   DTU Informatics, Technical University of Denmark. 

% Copyright Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

[I,J] = size(X);

opts = mgetopt(varargin);
alpha = updim(mgetopt(opts, 'alpha', zeros(I,N)),[I N]);
beta = updim(mgetopt(opts, 'beta', zeros(N,J)),[N J]);
theta = mgetopt(opts, 'theta', 0);
k = mgetopt(opts, 'k', 0);
A0 = mgetopt(opts, 'A', rand(I,N));
B0 = mgetopt(opts, 'B', rand(N,J));
sigma0 = mgetopt(opts, 'sigma', 1);
chains = mgetopt(opts, 'chains', 1);
stride = mgetopt(opts, 'stride', 1);
skip = mgetopt(opts, 'skip', 0);
nA = mgetopt(opts, 'nA', false(N,1));
nB = mgetopt(opts, 'nB', false(N,1));
ns = mgetopt(opts, 'ns', false);

Am = zeros(I,N,M*chains);
Bm = zeros(N,J,M*chains);
sigmam = zeros(M*chains,1);

x = sum(X(:).^2)/2;
for r = 1:chains
    A = A0;
    B = B0;
    sigma = sigma0;
    for m = 1:M
        for i = 1:skip*(m==1)+stride*(m>1)
            C = B*B';
            D = X*B';
            for n = 1:N
                if ~nA(n)
                    nn = [1:n-1 n+1:N];
                    A(:,n) = randr((D(:,n)-A(:,nn)*C(nn,n))/C(n,n), ...
                        sigma/C(n,n), alpha(:,n));
                end
            end
            if ~ns
                sigma = 1/gamrnd((I*J)/2+k, ...
                    1/(x+theta+sum(sum(A.*(A*C-2*D)))/2));
            end
            E = A'*A;
            F = A'*X;
            for n = 1:N
                if ~nB(n)
                    nn = [1:n-1 n+1:N];
                    B(n,:) = randr((F(n,:)-E(n,nn)*B(nn,:))'/E(n,n), ...
                        sigma/E(n,n), beta(n,:)');
                end
            end
        end
        Am(:,:,m+(r-1)*M) = A;
        Bm(:,:,m+(r-1)*M) = B;
        sigmam(m+(r-1)*M) = sigma;
    end
end

%--------------------------------------------------------------------------
function X = updim(x, dim)
% UPDIM Replicate and tile an array
%  
% Usage
%   X = updim(x, dim)

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

dimx = zeros(size(dim));
for k = 1:length(dim), dimx(k) = size(x,k); end
dim(dim==dimx) = 1;
X = repmat(x,dim);


%--------------------------------------------------------------------------
function x = randr(m, s, l)
% RANDR Random numbers from 
%   p(x)=K*exp(-(x-m)^2/s-l'x), x>=0 
%
% Usage
%   x = randr(m,s,l)

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

A = (l.*s-m)./(sqrt(2*s));
a = A>26;
x = zeros(size(m));

y = rand(size(m));
x(a) = -log(y(a))./((l(a).*s-m(a))./s);

R = erfc(abs(A(~a)));
x(~a) = erfcinv(y(~a).*R-(A(~a)<0).*(2*y(~a)+R-2)).*sqrt(2*s)+m(~a)-l(~a).*s;

x(isnan(x)) = 0;
x(x<0) = 0;
x(isinf(x)) = 0;
x = real(x);

%--------------------------------------------------------------------------
function out = mgetopt(varargin)
% MGETOPT Parser for optional arguments
% 
% Usage
%   Get alpha parameter structure from 'varargin'
%     opts = mgetopt(varargin);
%
%   Get and parse alpha parameter:
%     var = mgetopt(opts, varname, default);
%        opts:    parameter structure
%        varname: name of variable
%        default: default value if variable is not set
%
%     var = mgetopt(opts, varname, default, command, argument);
%        command, argument:
%          String in set:
%          'instrset', {'str1', 'str2', ... }
%
% Example
%    function y = myfun(x, varargin)
%    ...
%    opts = mgetopt(varargin);
%    parm1 = mgetopt(opts, 'parm1', 0)
%    ...

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

if nargin==1
    if isempty(varargin{1})
        out = struct;
    elseif isstruct(varargin{1})
        out = varargin{1}{:};
    elseif isstruct(varargin{1}{1})
        out = varargin{1}{1};
    else
        out = cell2struct(varargin{1}(2:2:end),varargin{1}(1:2:end),2);
    end
elseif nargin>=3
    opts = varargin{1};
    varname = varargin{2};
    default = varargin{3};
    validation = varargin(4:end);
    if isfield(opts, varname)
        out = opts.(varname);
    else
        out = default;
    end
    
    for narg = 1:2:length(validation)
        cmd = validation{narg};
        arg = validation{narg+1};
        switch cmd
            case 'instrset',
                if ~any(strcmp(arg, out))
                    fprintf(['Wrong argument %sigma = ''%sigma'' - ', ...
                        'Using default : %sigma = ''%sigma''\n'], ...
                        varname, out, varname, default);
                    out = default;
                end
            case 'dim'
                if ~all(size(out)==arg)
                    fprintf(['Wrong argument dimension: %sigma - ', ...
                        'Using default.\n'], ...
                        varname);
                    out = default;
                end
            otherwise,
                error('Wrong option: %sigma.', cmd);
        end
    end
end
