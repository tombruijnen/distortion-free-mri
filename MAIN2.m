%% DSPI method for distortion free MRI
clear all;close all;clc;
% Windows
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\MReconUMC\MReconUMC\Packages\utils'))
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\MReconUMC\MReconUMC\Standalone\Greengard_nufft\'))
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\Distortion_free_mri\distortion-free-mri\distortion-free-mri\'));
% cd('C:\Users\tombruijnen\Documents\Programming\MATLAB\RECON_MAIN\Main_Recon_Library\')

% Linux
addpath(genpath('/nfs/bsc01/researchData/USER/tbruijne/Projects_Main/Distortion_Free_MRI'))
cd('/nfs/rtsan02/userdata/home/tbruijne/Documents/Main_Recon_Library')
setup;

%% Body
% Set phantom and acquisition
ph.vis=1;
ph.woff_fat=220; % Hz
ph.N=256; % matrix size
ph.wf=.2;
ph.ff=.3;
ph.width_sus=4; % pixels
ph.B0=3; % Tesla
ph.gamma=2*pi*42.58e6;
ph.off_res_max=120; % Hz
ph.voxel_size=[1 1 1]; % mm
ph.R=100; % Undersampling factor R --> R^2 samples per time-point
ph.te=1.4E-03;
ph.adc_dur=2E-03;
ph.adc_dt=ph.adc_dur/ph.N;
ph.Id=[ph.N ph.N 1 1 ph.N];
ph.Kd=[ph.N round(ph.N/ph.R) 1 1 ph.N];

% Generate phantom + B0 maps
[F,W,B0,der_B0]=create_susceptibillity_phantom(ph);

% Generate k-space images and Cartesian distorted reconstruction
[cartI,Ist]=generate_cartesian_acquisition(F,W,B0,der_B0,ph);

% Create radial through time trajectory per time-point and density
[k,dcf]=radial_through_time_trajectory(ph);

% Create operators and store in struct op
op.W=DCF(dcf);
op.F=GG2D(k,ph.Kd);
op.S=CC(ones(ph.Id(1:4))); % no coil maps at the moment
op.A=PH(B0,ph.Id); % Known susceptibillity map, not true in vivo

% Retrieve radial k-space
rad_kspace=op.F*Ist;

% Check radial reconstruction
radI=op.F'*(op.W*(rad_kspace));

%% Low rank reconstruction
% SVD of true object
[U,S,V]=svd(reshape(Ist,[ph.N^2 ph.N]),'econ');S=diag(S);
figure,
subplot(4,4,[1 2 5 6]);plot(S,'LineWidth',2);title('SVD of object')
xlabel('Singular values');set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2)

% Compressed sensing recon 
par.TV=[0 0 0 0 0.1]; % lambdas in dimensions 
par.wavelet=0;
par.traj=k(:,:,:,:,:);
par.Niter=100;
par.kspace_data=single(rad_kspace(:,:,:,:,:));
par.csm=ones(ph.N,ph.N,'single');
CS=configure_compressed_sense(par,'bart');   

% Low rank reconstruction via BART
par=rmfield(par,'TV');
par.LR=[0 0 0 0 .2];
LR=configure_compressed_sense(par,'bart');   

close all;
figure,imshow3(squeeze(cat(5,demax(real(radI(:,:,:,:,1:50:end))),demax(real(CS(:,:,:,:,1:50:end))),demax(real(LR(:,:,:,:,1:50:end))),demax(real(Ist(:,:,:,:,1:50:end))))),[],[4 6])
figure,imshow3(squeeze(cat(5,demax(abs(radI(:,:,:,:,1:50:end))),demax(abs(CS(:,:,:,:,1:50:end))),demax(abs(LR(:,:,:,:,1:50:end))),demax(abs(Ist(:,:,:,:,1:50:end))))),[],[4 6])