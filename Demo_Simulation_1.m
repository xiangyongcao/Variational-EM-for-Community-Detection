% Demo of Simulation 1
clc;clear;

currentFolder = pwd;
addpath(genpath(currentFolder))

% the parameters of generating an attributed graph
Q = 2;
n = 100;
True_Alpha = [0.4,0.6];
True_PI=[0.8,0.05;0.05,0.4];
AtSize=[2,3];
True_Theta = cell(1,2);
True_Theta{1}=[0.7,0.3;0.8,0.2];
True_Theta{2}=[0.6,0.1,0.3;0.5,0.2,0.3];

time = 50;  % repeated number

% Store result for each time
PI11 = zeros(1,time);  PI12 = zeros(1,time); PI22 = zeros(1,time);
Alpha1 = zeros(1,time);
Theta11_1 = zeros(1,time); Theta21_1 = zeros(1,time);
Theta11_2 = zeros(1,time); Theta12_2 = zeros(1,time);
Theta21_2 = zeros(1,time); Theta22_2 = zeros(1,time);
CERR = zeros(1,time);
Modularity = zeros(1,time);
Entropy1 = zeros(1,time); Entropy2 = zeros(1,time);
Time = zeros(1,time);

for i = 1:time
    disp(['The ',num2str(i),'-th time is running.'])
    
    % Generate an attributed graph
    [X,Y,Glabel] = GGraph(n,True_Alpha,True_PI,True_Theta);
    % Call function
    [PI, Alpha, Theta, Tau, CluResult,modularity,entropy,time,ICL]=VEM(X, Y, Q, 50,'spectral');
    
    PI11(i) = max(diag(PI));PI12(i)=PI(1,2); PI22(i) = min(diag(PI));
    Alpha1(i) = max(Alpha);
    Theta11_1(i) = max(Theta{1}(1,:)); Theta21_1(i) = max(Theta{1}(2,:));
    Theta11_2(i) = max(Theta{2}(1,:)); Theta12_2(i) = min(Theta{2}(1,:));
    Theta21_2(i) = max(Theta{2}(2,:)); Theta22_2(i) = min(Theta{2}(2,:));
    CERR(i) = CER(Glabel,CluResult);
    Modularity(i) = modularity;
    Entropy1(i) = entropy(1); Entropy2(i) = entropy(2);
    Time(i) = time;
end
disp(['The program stops!'])
% Statistics of the results
mean_PI11 = mean(PI11); std_PI11 = std(PI11);
mean_PI12 = mean(PI12); std_PI12 = std(PI12);
mean_PI22 = mean(PI22); std_PI22 = std(PI22);
mean_Alpha1 = mean(Alpha1); std_Alpha1 = std(Alpha1);
mean_Theta11_1 = mean(Theta11_1);  std_Theta11_1 = std(Theta11_1);
mean_Theta21_1 = mean(Theta21_1); std_Theta21_1 = std(Theta21_1);
mean_Theta11_2 = mean(Theta11_2); std_Theta11_2 = std(Theta11_2);
mean_Theta12_2 = mean(Theta12_2); std_Theta12_2 = std(Theta12_2);
mean_Theta21_2 = mean(Theta21_2); std_Theta21_2 = std(Theta21_2);
mean_Theta22_2 = mean(Theta22_2); std_Theta22_2 = std(Theta22_2);
mean_CERR = mean(CERR); std_CERR = std(CERR);
mean_Modularity = mean(Modularity); std_Modularity = std(Modularity);
mean_Entropy1 = mean(Entropy1); std_Entropy1 = std(Entropy1);
mean_Entropy2 = mean(Entropy2); std_Entropy2 = std(Entropy2);
mean_Time = mean(Time); std_Time = std(Time);

PI11 = [mean_PI11,std_PI11];
PI12 = [mean_PI12,std_PI12];
PI22 = [mean_PI22,std_PI22];
Est_PI = [mean_PI11,mean_PI12;mean_PI12,mean_PI22];
Alpha1 = [mean_Alpha1, std_Alpha1];
Theta11_1 = [mean_Theta11_1,std_Theta11_1];
Theta21_1 = [mean_Theta21_1,std_Theta21_1];
Theta11_2 = [mean_Theta11_2,std_Theta11_2];
Theta12_2 = [mean_Theta12_2,std_Theta12_2];
Theta21_2 = [mean_Theta21_2,std_Theta21_2];
Theta22_2 = [mean_Theta22_2,std_Theta22_2];
CER = [mean_CERR,std_CERR]; 
modularity = [mean_Modularity,std_Modularity]; 
entropy = [mean_Entropy1,std_Entropy1;mean_Entropy2,std_Entropy2];
time = [mean_Time std_Time];
disp(['True Alpha1 is ']);
True_Alpha1 = max(True_Alpha);
True_Alpha1

disp(['Estimated Alpha1 is ']);
mean_Alpha1

disp(['True PI is ']);
True_PI

disp(['Estimated PI is ']);
Est_PI

