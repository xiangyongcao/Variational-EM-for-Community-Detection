% Generate an attributed graph

function [X,Y,Glabel] = GGraph(n,Alpha,PI,Theta)

% Q = length(Alpha);
P = length(Theta);

% Generate label for each vertex
Z = mnrnd(1,Alpha,n);

Glabel = zeros(1,n);
for i = 1:n
    Glabel(i) = find(Z(i,:)==1);
end

% Generate adjacent matrix
B = Z * PI * Z';
C = rand(n);
X = (C<=B);
X = X - diag(diag(X));

% Generate attributed matrix
Y = zeros(n,P);
for p = 1:P
   theta = Theta{p};
   for i = 1:n
       alpha = theta(Glabel(i),:);
       z = mnrnd(1,alpha,1);
       Y(i,p) = find(z==1);
   end
end



% for i = 1:n
%     for j = 1:n
%       edge_seed = rand(1);
%       if P(Glabel(i),Glabel(j)) <= edge_seed
%           A(i,j) = 1;
%       end
%     end
% end


