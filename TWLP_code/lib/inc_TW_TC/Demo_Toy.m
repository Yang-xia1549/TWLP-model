function [MatrixTimes,MatrixResults,Re_tensor,Out,G,C,B]=Demo_Toy(image,SR, R,rho,alpha,miu,omiga,p,lamda,rho_rate,eta,Omega)

% clc;
% clear;
% close all;
addpath(genpath([pwd,'/lib']));
addpath(genpath([pwd,'/data']));
%%
%load('CAVE_Toy.mat')
X=image;
% X(X<0) = 0;
% if max(X(:))>1
%     X = X/max(X(:));%归一化
% end
%% Generage corrupted image
% maxP = max(abs(X(:)));
% [n1,n2,n3] = size(X);
% Xn = X;
% rhos = 0.1;                       % setting corrupted level eta=10 
% rng(234)
% ind = find(rand(n1*n2*n3,1)<rhos);
% Xn(ind) = rand(length(ind),1);    % noisy samples      
% [n1,n2,n3] = size(X);    % sample size:481X321X3 or 321X481X3
% Xn=Noisy(X,'Random',noise);
%% 采样
sample_ratio = SR;
fprintf('### Performing SR: %4.2f ###\n', sample_ratio);
Y_tensor = X;
clear X
Ndim = ndims(Y_tensor);
Nway = size(Y_tensor);
% rng(2)
% Omega = find(rand(prod(Nway),1)<sample_ratio);
Y_init = initialization_M(Nway, Omega, Y_tensor(Omega));
%%
Algorithms = {'Observed','TW-TC'};
EN_TW_TC = 1;
A_num = length(Algorithms);
Re_tensor = cell(A_num,1);
NumIndexes = 2;  
MatrixResults = zeros(A_num, NumIndexes);
MatrixTimes = zeros(A_num, 1);
alg = 0;
fprintf('###################### Please wait......######################\n')
%% Miss_Data
alg = alg+1;
F = zeros(Nway);
F(Omega) = Y_tensor(Omega);
Re_tensor{alg} = F;
[MatrixResults(alg,:)] = quality_access(Re_tensor{alg}*255,image*255);
%% Perform Our TW-TC Algorithm
alg = alg+1;
% start=0.1;
% end_value=3;
% step=0.1;
% i=1;
% for omiga=start:step:end_value
%     i=i+1;
if EN_TW_TC
     fprintf('\n');
    disp(['performing ', Algorithms{alg}, ' p', num2str(p), ' omiga', num2str(omiga), ' lamda ', num2str(lamda),' miu',num2str(miu),' eta',num2str(eta)]);
    opts = [];
    opts.tol   = 1e-5;
    opts.maxit = 200;
    opts.rho   = rho;%1.9;%1.01
    opts.Xtrue = image;
    opts.R =R;% [5, 17, 5; % R_i
             % 3, 3, 3]; % L_i
    opts.alpha = alpha;%[ 1e-4 1.7 1e-4 1e-4];
    opts.miu=miu;
    opts.omiga=omiga;
    opts.lamda=lamda;
    opts.eta=eta;
    opts.p=p;
    opts.rho_rate=rho_rate;
    t_s = tic;
    [Re_tensor{alg}, G, C, Out,B] = inc_TWLP_TC(Y_init, Omega, opts);
    MatrixTimes(alg,:) = toc(t_s);
    [MatrixResults(alg,:)] = quality_access(Re_tensor{alg}*255,image*255);
end
% end
%
fprintf('\n');
    disp([MatrixResults(2,1)]);
fprintf('###################### Complete execution! ! !######################\n')
disp([' MatrixResults', num2str(MatrixResults(2,1))]);
 fprintf('\n');
%% Show result
% location  = [10 50 45 85];
% range_bar = [0, 0.2];
% save_dir = 'D:\桌面\111\MSI -noisy\lib\result\NEW\';
% %
% filename1=['image_' image_name '_SR_ ' num2str(SR) ' _R1_' num2str(R(1,1)) '_R2_' num2str(R(1,2)) '_R3_' num2str(R(1,3))   '_residual.jpg'];
% filename=['image_' image_name '_SR_ ' num2str(SR) ' _R1_' num2str(R(1,1)) '_R2_' num2str(R(1,2)) '_R3_' num2str(R(1,3)) '_image.jpg'];
% Image_truncate(image, Re_tensor, location);    % Display the RGB images
% saveas(gcf, fullfile(save_dir, filename)); % 保存为JPG
% close(gcf);
% %
% Image_residual(image, Re_tensor, range_bar);   % Display the error maps 
% saveas(gcf, fullfile(save_dir, filename1));
% close(gcf);

end