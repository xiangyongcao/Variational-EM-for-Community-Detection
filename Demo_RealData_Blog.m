% Demo of real data: Political Blog
clc;
clear;
% load data
load('X_blog_directed.mat');
load('Y_blog.mat');
Q = 2;

% Call function
[PI, Alpha, Theta, Tau, CluResult,modularity,entropy,time,ICL]=VEM(X_blog,Y_blog,Q,50,'spectral','directed');

% Performance Assesment
cer = CER(Y_blog,CluResult);
disp(['modularity:',num2str(modularity)]);         
disp(['entropy:',num2str(entropy)]);           
disp(['CER:',num2str(cer)]);     
disp(['optimizing time(s):',num2str(time)]);