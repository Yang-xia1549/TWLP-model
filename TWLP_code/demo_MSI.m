clc;
clear;
close all;
addpath(genpath([pwd,'/lib']));
addpath(genpath([pwd,'/data']));


for SR=[0.1]

    image_name='beens';

    load([image_name '.mat']);
    image=double(image);


    Nway=size(image);
    D_dim=ndims(image);
    rng(2)
    Omega = find(rand(prod(Nway),1)<SR);
    Y_init = initialization_M(Nway, Omega, image(Omega));

    p=0.6;
    F = zeros(Nway);
    F(Omega) = image(Omega);
    Re_tensor{1} = F;
    [MatrixResults(1,:)] = quality_access(Re_tensor{1}*255,image*255);
    %% TWLP
    alpha=zeros(D_dim,4);
    alpha_t=[1 1 0 0];
    for L1=3
        for R1=2
            for R2=25
                for R3=5
    for omiga=1 
        for miu=0.01
            for eta=10 
                for lamda=2 
                    for rho=0.1 %[0.01 0.01 0.1]
                        R=[R1, R2, R3;
                            L1, L1,R3];
                        for k=1:D_dim
                            alpha(k,:)=[R(2,k) R(1,k) R(2,k) R(1,k)];
                            alpha(k,:)=alpha(k,:).*alpha_t;
                        end
                        fprintf('\n');
                        disp(['performing p', num2str(p), ' omiga', num2str(omiga), ' lamda ', num2str(lamda),' miu',num2str(miu),' eta',num2str(eta)]);
                        opts = [];
                        opts.tol   = 1e-5;
                        opts.maxit = 200;
                        opts.rho   = rho;
                        opts.Xtrue = image;
                        opts.R =R;
                        opts.alpha = alpha;
                        opts.miu=miu;
                        opts.omiga=omiga;
                        opts.lamda=lamda;
                        opts.eta=eta;
                        opts.p=p;
                        t_s = tic;
                        [Re_tensor{2}, G, C, Out,B] = inc_TWLP_TC(Y_init, Omega, opts);
                        MatrixTimes(2,:) = toc(t_s);
                        [MatrixResults(2,:)] = quality_access(Re_tensor{2}*255,image*255);
                        rse = perfscore(Re_tensor{2}*255,image*255);
                        fprintf('\n');
                        disp([MatrixResults(2,1)]);
                        fprintf('###################### Complete execution! ! !######################\n')
                    end
                end
            end
        end
    end
                end
            end
        end
    end

    %% Show result
location  = [10 50 45 85];
range_bar = [0, 0.2];
Image_truncate(image, Re_tensor, location);    % Display the RGB images
Image_residual(image, Re_tensor, range_bar);   % Display the error maps 
end
