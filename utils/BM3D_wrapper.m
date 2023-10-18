function [image_fasthyde] = BM3D_wrapper(image_ori,SIGMA)
[Lines, Columns, B] = size(image_ori);
N=Lines*Columns;
for i=1:size(image_ori,3)
    Y(i,1:N)= reshape(image_ori(:,:,i),[1,N]);
end

%% denoising:

y_est_bm3d = zeros(Lines,Columns,B);

eigen_Y_bm3d=[];

for i=1:B
    % produce eigen-image
    eigen_im = Y(i,:);
    min_x = min(eigen_im);
    max_x = max(eigen_im);
    eigen_im = eigen_im - min_x;
    scale = max_x-min_x;
    
    %scale to [0,1]
    eigen_im = reshape(eigen_im, Lines, Columns)/scale;
    sigma =  SIGMA(i)/scale;
    filt_eigen_im = eigen_im;
    
    if sigma>0
        
        [dummy, filt_eigen_im] = BM3D(1,eigen_im, sigma*255);
    end
    
    eigen_Y_bm3d(i,:) = reshape(filt_eigen_im*scale + min_x, 1,N);
end

% reconstruct data using denoising engin images
Y_reconst_bm3d = eigen_Y_bm3d;

image_fasthyde=[];
for i=1:B
    image_fasthyde(1:Lines,1:Columns,i) = reshape(Y_reconst_bm3d(i,:),[Lines,Columns]);
end

t2=clock;

end