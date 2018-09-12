%% Investigate realistic low rankness of the off-resonance dimension
clear all;close all;clc;
% Windows
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\MReconUMC\MReconUMC\Packages\utils'))
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\MReconUMC\MReconUMC\Standalone\Greengard_nufft\'))
% addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\Distortion_free_mri\distortion-free-mri\distortion-free-mri\'));
% cd('C:\Users\tombruijnen\Documents\Programming\MATLAB\RECON_MAIN\Main_Recon_Library\')

% Linux
addpath(genpath('/nfs/bsc01/researchData/USER/tbruijne/Projects_Main/Distortion_Free_MRI'))
addpath(genpath('/nfs/bsc01/researchData/USER/tbruijne/Projects_Software/XCAT/'))
cd('/nfs/rtsan02/userdata/home/tbruijne/Documents/Main_Recon_Library')
setup;

%% General settings
ph.R=100; % Undersampling factor R --> R^2 samples per time-point
ph.N=256; % matrix size (2D)
ph.voxel_size=[1.5 1.5 2]; % [mm]
ph.matrix=[ph.N ph.N 1 1 ph.N];
ph.B0=3; % 3T
ph.woff_fat=73*ph.B0;  % Off-resonance of fat
ph.ff=.2; % fat fraction
ph.gamma=2*pi*42.58e6;
ph.width_sus=15; % Width of titanium alloy [pixels]
ph.off_res_max=500; % Maximum off-resonance of linear ramped field [Hz]
ph.te=1.4E-03;
ph.adc_dur=2E-03;
ph.adc_dt=ph.adc_dur/ph.N;
ph.Id=[ph.N ph.N 1 1 ph.N];
ph.Kd=[ph.N round(ph.N/ph.R) 1 1 ph.N];
ph.noise_levels=10;

%% Extract abdominal, lung and brain slice and investigate low-rankness on multiple noise levels. 
xcat_path='/nfs/bsc01/researchData/USER/tbruijne/Projects_Software/XCAT/XCAT_V2_LINUX';

% Generate abdomen\lung\brain
clear pars
pars.pixel_width=10^-1*ph.voxel_size(1);   % X/Y resolution [cm]
pars.slice_width=10^-1*ph.voxel_size(3);   % Z resolution [cm]
pars.array_size=ph.N;   % X/Y matrix size [-]
pars.startslice=600;     % Z cropping begin [-]
pars.endslice=865;     % Z cropping end [-]
xcat=XCAT(xcat_path,[],pars,'MR');
abdo=xcat(:,:,1);
lung=xcat(:,:,100);
brain=xcat(:,:,255);clear xcat

%% Add off-resonance map, water/fat compartment and titanium insert
% Water/Fat separation
[Wbrain,Fbrain]=multi_compartmentalise(brain,ph.ff);
[Wabdo,Fabdo]=multi_compartmentalise(abdo,ph.ff);
[Wlung,Flung]=multi_compartmentalise(lung,ph.ff);

% Generate off-resonance fields with titanium inserts
[Fbrain,Wbrain,B0,der_B0]=create_susceptibillity_phantom(Wbrain,Fbrain,ph);
[Fabdo,Wabdo,B0,der_B0]=create_susceptibillity_phantom(Wabdo,Fabdo,ph);
[Flung,Wlung,B0,der_B0]=create_susceptibillity_phantom(Wlung,Flung,ph);

% Generate temporal k-space and Cartesian distorted reconstruction
[~,Ibrain]=generate_cartesian_acquisition(Wbrain,Fbrain,B0,der_B0,ph);
[~,Iabdo]=generate_cartesian_acquisition(Wabdo,Fabdo,B0,der_B0,ph);
[~,Ilung]=generate_cartesian_acquisition(Wlung,Flung,B0,der_B0,ph);

%% SVD analysis for the multiple anatomical sites and different SNR levels
noise_sigma=linspace(0,40,ph.noise_levels);
noise_sigma = 0;
for n=1:ph.noise_levels
    fnoise=@(x)(x+(-.5*noise_sigma(n)+(noise_sigma(n)*rand(size(x))))+1j*(-noise_sigma(n)/2+noise_sigma(n)*rand(size(x))));    
    [U,S,V]=svd(reshape(fnoise(Ibrain),[ph.N^2 ph.N]),'econ');Sbrain(:,n)=diag(S);
    [U,S,V]=svd(reshape(fnoise(Iabdo),[ph.N^2 ph.N]),'econ');Sabdo(:,n)=diag(S);
    [U,S,V]=svd(reshape(fnoise(Ilung),[ph.N^2 ph.N]),'econ');Slung(:,n)=diag(S);
    n
end
figure,
subplot(221);
imshow3(cat(3,brain,abdo,lung),[],[1 3])
subplot(222);
plot(Sbrain,'k','LineWidth',2);hold on;plot(Sabdo,'b','LineWidth',2);plot(Slung,'r','LineWidth',2);
xlabel('Singular values');legend('Brain','Abdo','Lung');title('SVD of temporal image');axis([0 30 0 10*10^4])
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2);grid on