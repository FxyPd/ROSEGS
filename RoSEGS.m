function [EZ_est] =  RoSEGS(img_noisy,k_subspace,img_clean,alpha,beta,rho)
% addpath(genpath('E:\Fuxiyou\Hyperspectral\super-resolution\CNN-FUS-master\CNN_subspace\FFDNet-master'));
%Input:
% img_noisy          noisy 3D data of size row*column*band
% k_subspace         dimension of subspace.
% img_clean          clean 3D data of size row*column*band.
% alpha              alpha
% beta               beta
% rho                rho

% Output:
% EZ_est denoised 3D image


[row, column, band] = size(img_noisy);
[nr, nc, band] = size(img_noisy);

N=row*column;
Y_noisy = reshape(img_noisy, [N band])';

%  -------------Subspace Learning Against Mixed Noise---------------------
%An adaptive median filter is applied to noisy image to remove the bulk of
%impulse noise and stripes
if 1
    for in=1:N
        % hampel filtering
        [Y_median(:,in),Y_position(:,in),xmedian,xsigma] = hampel(Y_noisy(:,in),3,3);
    end
    img_median = reshape(Y_median', row, column, band);
end
img_remove_outlier = img_median;
Y_remove_outlier = reshape(img_remove_outlier, N, band)';
[~,Rw] = estNoise(Y_remove_outlier,'additive');
Rw_ori = Rw;

%data whitening so that noise variances of each band are same
Y_noisy = inv(sqrt(Rw))*Y_noisy;
Y_median = inv(sqrt(Rw))*Y_median;
Y_remove_outlier= inv(sqrt(Rw))*Y_remove_outlier;

% simulate the removal of some bad pixels
% ind=randperm(size(Y_remove_outlier,2));
% percent=round(size(Y_remove_outlier,2)*0.0000001);% 0.1 is pecentage of bad pixels
% Y_remove_outlier(:,ind(1:percent))=[];

%Subspace learning from the coarse image without stripes and impulse noise
[E,S,~] = svd(Y_remove_outlier*Y_remove_outlier'/N);
E = E(:,1:k_subspace);
p = k_subspace;


Z = E'*Y_median;
%% initializtion:
mu=10;
mu1=1e-2;mu2=1;mu3=1e-2;mu4=1e-2;% 


V1=zeros(size(Y_noisy));
G1=zeros(size(V1));

if 1
M=double(~Y_position);
nn=size(M,1)*size(M,2);
% ee=ones(nn,1)/1.1;
ee=1./(ones(nn,1)+reshape(M,nn,1)*1/mu1);
diagMM = spdiags([ee-ee ee ee-ee],-1:1,nn,nn);
end


V2=Z;
G2=zeros(size(Z));

Eny_x_fft   = (abs(psf2otf([+1; -1], [nr,nc,p]))).^2  ;
Eny_y_fft   = (abs(psf2otf([+1, -1], [nr,nc,p]))).^2  ;
Eny_fft  =  Eny_x_fft + Eny_y_fft;

% auxiliary functins
Vec = @(x)(x(:));
inVec = @(x)(reshape(x,size(Y_noisy,1),size(Y_noisy,2)));

Rx = zeros([nr,nc,p]);      % R1 : auxiliary variable for GS_x
Ry = zeros([nr,nc,p]);      % R2 : auxiliary variable for GS_y
G3x = Rx;                % W3 : multiplier for DQ_x-R1
G3y = Ry;                % W4 : multiplier for DQ_y-R2

Q = Z;                  % Q :  auxiliary variable for A
Q=hyperConvert3D(Q,nr, nc );
G4=zeros([nr,nc,p]);

k=0;

Y=Y_noisy;
%% iteration
for iter = 1: 30
    
    k = k+1
    %% update V1
    W=(Y-E*Z);
    
    V1=diagMM*Vec(W);

    V1=inVec(V1);
    
    
    %% update A
    tmpZ=mu1*(E'*E)+mu2*eye(size(E'*E,1))+mu4*eye(size(E'*E,1));
    Z=tmpZ\(mu1*E'*(Y-V1-G1/2/mu1)+mu2*(V2+G2/2/mu2)+mu4*hyperConvert2D(Q+G4/2/mu4));
    
    Zt=hyperConvert3D(sqrt(Rw_ori)*E*Z,nr,nc);
    
    psnr = PSNR_m(img_clean,Zt);
    psnr.ave
       
    %% update V2 BM3D denoiser
    if 1 %BM3D
        B2=Z-G2/(2*mu2);
        B2=hyperConvert3D(B2,nr, nc );
        noise_std = ones(p,1)*alpha/2/mu2;%alpha/2/mu2
        [V2] = BM3D_wrapper(B2,noise_std);
        V2=hyperConvert2D(V2);
    end
    
    %% - Q subproblem update
    diffT_p = diffT3(2*mu3 * Rx - G3x, 2*mu3 * Ry - G3y, [nr,nc,p]);%W4 or W3?
    temp1 = reshape(diffT_p + 2*mu4*hyperConvert3D(Z,nr,nc) - G4,  [nr,nc,p]);
    z = real(ifftn(fftn(temp1) ./ (2*mu3*Eny_fft + 2*mu4)));
    Q = reshape(z, [nr,nc,p]);
    
    %% - R1 and R2 subproblem update
    
    [diff_Qx, diff_Qy] = diff3(Q,  [nr,nc,p]);
    Rx = Thres_21(diff_Qx+ G3x/mu3/2, beta/mu3/2);%beta/mu3/2 0.06
    Ry = Thres_21(diff_Qy+ G3y/mu3/2, beta/mu3/2);%beta/mu3/2 0.06
    
    %% update E
    E_est   = -(V1-Y+G1/2/mu1)*Z';
    [U,~,V] = svd(E_est,'econ');
    E = U*V';
    
    %% updata Lagarigian parameters
    G1  = G1 +2*mu1*(V1-Y+E*Z);
    G2  = G2 +2*mu2*(V2-Z);
    G3x = G3x+2*mu3*(diff_Qx-Rx);
    G3y = G3y+2*mu3*(diff_Qy-Ry);
    G4=G4+2*mu4*(Q-hyperConvert3D(Z,nr,nc));
    

    mu2=min(mu,mu2*rho);

end

EZ_est = hyperConvert3D(sqrt(Rw_ori)*E*Z,nr,nc);


