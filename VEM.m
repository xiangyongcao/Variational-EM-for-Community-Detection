function [PI, Alpha, Theta, Tau, CluResult,modularity,entropy,time,ICL] = ...
    VEM(X, Y, Q, MaxIter,IniType, NetType)
%
% Author - Xiangyong Cao, 05/2014
%
% Email  - caoxiangyong45@gmail.com
%
% Description - main function to implement the variational EM algorithm.
%
% Input  - X    : adjacency matrix   n x n
%        - Y    : attribute matrix   n x P
%        - Q    : cluster number
%        - MaxIter   : iteration number, default 50
%        - IniType   : 'random','spectral','Int-Cluster'
%        - NetType   : 'undirected','directed'
% Output - PI        : connection matrix   Q x Q
%        - Alpha     : proportion of each group
%        - Theta     :
%        - Tau       : posterior probability of each vertex
%        - CluResult : label of each vertex
%        - modularity : scalar, structure quality measure
%        - entropy    : Tx1 entropy vector with an element for each attribute;
%                       attribute quality measure.
%        - time       : optimizing time
% -------------------------------------------------------------------------

% default settings-----------------------------


if nargin < 6
    NetType = 'undirected';
end
if nargin < 5
    IniType = 'Inc-Cluster';
end
if nargin < 4
    MaxIter = 50;
end

% common constants---------------------
[n,P] = size(Y);
AtSize = zeros(1,P);
for p =1:P
    AtSize(p) = length(unique(Y(:,p)));
end


% format data --------------------------------------
Y1 = formatted(n,Y);

% Init --------------------------------------------
switch IniType
    case 'random'
        tmp = rand(n,Q);
        Tau = tmp./repmat(sum(tmp,2),1,Q);
    case 'spectral'
        [ C, L, U ] = SpectralClustering( X, Q, 3 );
        Tau = full(C);
    case 'Inc-Cluster'
        posterior = init(X,Y1,Q);
        Tau = full(posterior);
end

% ----------internal control parameters----------
int_maxIter = 5;
out_tol = 1e-10;
int_tol = 1e-4;
log_min = log(realmin);
% prune   = 1e-10;

tic;
out_iter = 1;
Theta = cell(P,1);
while out_iter<MaxIter
    
    out_iter;
    
    %% Variational M step
    % update Alpha, PI, Theta
    col_sum = sum(Tau);
    Alpha = col_sum/n;
    PI = (Tau' * X * Tau)./(col_sum'*col_sum - Tau'*Tau);
    
    for p=1:P
        Theta{p} = Tau'*Y1{p};
        Theta{p} = Theta{p}./repmat(sum(Theta{p},2),1,AtSize(p));
    end
    
    % convergence
    if out_iter > 1
        old_L = L;
    end
    L = lblikelihood(X,Y1,Tau,Alpha,PI,Theta,NetType);
    if out_iter > 1
        err = abs(1-L/old_L);  % relative error of lower bound
        if err < out_tol
            break;
        end
    end
    
    
    %% Variational E step
    % updata Tau
    A = repmat(log(Alpha),n,1);
    for p = 1:P
        a = log(Theta{p});
        A = A + Y1{p}*a';
    end
    
    b = log(1-PI);
    B1 = log(PI) - b;
    int_iter = 1;
    diff = int_tol + 1;
    keep_on = true;
    while (int_iter<int_maxIter) && (diff>int_tol) &&(keep_on)
        old_Tau = Tau;
        tmp = Tau*b';
        tmp = repmat(sum(tmp,1),n,1)-tmp;
        tmp1 = X*Tau*B1';
        switch NetType
            case 'undirected'
                Tau = 0.5*(tmp1+tmp)+A;
            case 'directed'
                Tau= tmp1+tmp+A;
        end
        Tau(Tau<log_min) = log_min;
        Tau = Tau-repmat(max(Tau,[],2),1,Q);
        Tau = exp(Tau);
        Tau = Tau./repmat(sum(Tau,2),1,Q);
        %         Tau(Tau<prune) = 0;
        if int_iter > 1
            old_diff = diff;
            diff = max(max(abs(Tau-old_Tau)));
            if diff >= old_diff
                keep_on = false;    % diff increase, stop iterating
            end
        elseif int_iter == 1
            diff = max(max(abs(Tau-old_Tau)));
        end
        int_iter = int_iter + 1;
    end
    out_iter = out_iter + 1;
end
time = toc;

%CluResult
[~,CluResult1] = max(Tau,[],2);

if Q == 1
    CluResult = ones(1,n);
else
    CluResult = CluResult1;
end
[modularity,entropy] = evaluate(X,Y1,CluResult);

ICL = L - 0.5*(Q-1)*log(n)- 0.25*Q*(Q+1)*log(0.5*n*(n-1))-0.5*Q*log(n)*(sum(AtSize)-P);

