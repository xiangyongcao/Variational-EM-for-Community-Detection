 % compute the lower bound
function L = lblikelihood(X,Y1,Tau,Alpha,PI,Theta,NetType)

n = size(X,1);
P = length(Theta); 

tmp = sum(sum(Tau(Tau>0).*log(Tau(Tau>0))));
tmp1 = Tau.*repmat(log(Alpha),n,1);
tmp1(isnan(tmp1)) = 0;
C = sum(sum(tmp1)) - tmp;

B = 0;
for p = 1:P
   tmp = Tau.*(Y1{p}*log(Theta{p}'));
   tmp(isnan(tmp)) = 0;
   B = B + sum(sum(tmp));
end
b = log(1-PI);
a = log(PI) - b;
tmp = Tau*b';
tmp = repmat(sum(tmp),n,1) - tmp;
tmp1 = Tau.*(X*Tau*a');
tmp1(isnan(tmp1))=0;
tmp2 = Tau.*tmp;
tmp2(isnan(tmp2))=0;
A = sum(sum(tmp1))+sum(sum(tmp2));

switch NetType
    case 'undirected'
        L = 0.5*A + B + C;
    case 'directed'
        L = A + B + C;
end


