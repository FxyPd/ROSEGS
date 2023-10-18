%% Name: RoSEGS
%
%  Generate the denoising results of RoSEGS, please cite our paper when using these codes:
%
%  [1] X. Fu, Y. Guo, M. Xu and S. Jia, "Hyperspectral Image Denoising via 
%      Robust Subspace Estimation and Group Sparsity Constraint," in IEEE 
%      Transactions on Geoscience and Remote Sensing, vol. 61, pp. 1-16, 
%      2023, Art no. 5512716, doi: 10.1109/TGRS.2023.3277832.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT NOTE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      This script uses the package BM3D  (v2 (30 January 2014))
%      to implement the denoising algorithm BM3D introduced in
%
%      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image denoising by
%      sparse 3D transform-domain collaborative filtering," IEEE Trans.
%      Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%
%      The BM3D package  is available at the
%      BM3D web page:  http://www.cs.tut.fi/~foi/GCF-BM3D
%
%      Download this package and install it is the folder /BM3D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Xiyou Fu (fuxiyou@qq.com)
%         June, 2023
%%


clear;clc;close all;
% Q: How to run the different set of parameters?
% A: dataset              'Pavia' - Pavia image;
%                         'DC' - DC image.
%    case_num             '1' - non i.i.d. Guassian noise;
%                         '2' - non i.i.d. Guassian noise + stripes;
%                         '3' - non i.i.d. Guassian noise + salt&pepper noise;
%                         '4' - non i.i.d. Guassian noise + stripes + salt&pepper noise
%
% USAGE EXAMPLES:
% case_num = 4; %1:4
% dataset ='Pavia' ;
% %--------------------------
% case_num = 1;
% dataset ='DC';
% %--------------------------

addpath(genpath('/BM3D'));
addpath('./utils');
addpath('./data');

for case_num= 1%1:4
    dataset ='Pavia';
    %  dataset ='DC';
    switch case_num
        case 1  % non i.i.d. Guassian noise
            k_subspace = 8;
            k = 233;
            stripes = 0;
            impulse = 0;
            stripe_band_num = 30; % the number of band to add stripes. DC -> 60       Pavia -> 30
            impluse_ratio = 0.005;
            noise_simulation_Gaussian;   % use the script to simulate the synthetic noise
        case 2  % non i.i.d. Guassian noise + stripes
            k_subspace = 8;
            k = 233;
            stripes = 1;
            impulse = 0;
            stripe_band_num = 30; % the number of band to add stripes. DC -> 60       Pavia -> 30
            impluse_ratio = 0.005;
            noise_simulation_Gaussian;   % use the script to simulate the synthetic noise
        case 3  % non i.i.d. Guassian noise + salt&pepper noise
            k_subspace = 8;
            k = 233;
            stripes = 0;
            impulse = 1;
            stripe_band_num = 30; % the number of band to add stripes. DC -> 60       Pavia -> 30
            impluse_ratio = 0.005;
            noise_simulation_Gaussian;   % use the script to simulate the synthetic noise
        case 4  % non i.i.d. Guassian noise + stripes + salt&pepper noise
            k_subspace = 8;
            k = 233;
            stripes = 1;
            impulse = 1;
            stripe_band_num = 30; % the number of band to add stripes. DC -> 60       Pavia -> 30
            impluse_ratio = 0.005;
            noise_simulation_Gaussian;   % use the script to simulate the synthetic noise 
    end
    
    
    %% compute the quantitive assement indexes of the noisy HSI
    num = 1;
    methodname{num}= 'Noisy';
    disp('*********************** noisy ************************');
    [MPSNR(num),PSNRV(:,num),MSSIM(num),SSIMV(:,num),MFSIM(num),FSIMV(:,num) ] = QuanAsse_psnr_ssim_fsim(img_clean,img_noisy);

    %------------- Proposed -------------
    num=2;
    methodname{num}= 'Proposed';

    alpha=0.6;
    beta=2; % Note: the abscissa value in Fig.9 of RoSEGS paper is not correct. 
            % The correct value should be the abscissa value devided by 50, 
            % i.e., the optimal value of beta should be 2.    
    rho=1.04;% rho 
    
    t1=tic;
    [Y_proposed]=  RoSEGS(img_noisy,k_subspace,img_clean,alpha,beta,rho);%pa is set to 3;
    runingtime(num) = toc(t1);
    
    [MPSNR(num),PSNRV(:,num),MSSIM(num),SSIMV(:,num),MFSIM(num),FSIMV(:,num) ] = QuanAsse_psnr_ssim_fsim(img_clean,Y_proposed);
    
%     clear;clc;close all;
end

