%% Design through time radial k-space trajectory and corresponding trajectory
% |o----|  |-----|       |----o|    CENTER OUT RADIAL TRAJECTORY
% |-o---|  |-----|       |---o-|        - N time-points
% |--o--|  |--o--|  ...  |--o--|        - M reconstruction matrix
% |-----|  |-o---|       |-----|        - A radial arms
% |-----|  |o----|       |-----|        - Phi radial angle between arms
%
% Each radial arm has N samples per definition.
% Must have an integer number A 
% 
%
% T. Bruijnen @ 20180408
clear all;close all;clc
addpath(genpath('C:\Users\tombruijnen\Documents\Programming\MATLAB\Distortion_free_mri\distortion-free-mri\distortion-free-mri'))

%% Settings
par.fov=40;
par.N=15;
par.M=128; 
par.A=1;
par.maxk=1/par.fov*par.M;
par.Smax=20;
par.Gmax=2.5;
par.dt=1e-2;
par.ds=0.001;

% Default values
par.ga1=2*pi/par.N;
par.ga2=2*pi/par.A;
par.cco=linspace(.7,0,par.N)'*[1 1 1];

%% Calculate center-out radial k
par.ktraj=dspi_center_out_radial_trajectory(par);

%% Calculate time-optimal gradient waveforms
% Demodulate with k0
[time,wf]=calculate_demodulated_trajectory(par,1);









%% Get time-optimal gradient waveform designs

time=[];
% tic
% for ro=1:N
%     C=[real(ktraj(1:end,ro,1)) imag(ktraj(1:end,ro,1)) zeros(size(ktraj(1:end,ro,1)))];
%     [~,time_riv,g_riv(:,:,ro)]=minTimeGradient(C,0,g0,0,Gmax,Smax,dt);          % Rotationally invariant solution
%     time(:,ro)=linspace(0,time_riv,size(g_riv,1));
% end
% toc

C=[real(ktraj(1:end,nro,1)) imag(ktraj(1:end,nro,1)) zeros(size(ktraj(1:end,nro,1)))];
[~,time_riv,g_riv]=minTimeGradient(C,0,0,0,Gmax,Smax,dt,0.0001);          % Rotationally invariant solution
time=linspace(0,time_riv,size(g_riv,1));

subplot(223);plot(time,g_riv,'LineWidth',2);
axis([0 time(end) -1.5*Gmax 1.5*Gmax]);xlabel('Time [ms]');ylabel('G [cT/m]');
title('Example of k-space at t=tau');grid on;box on
set(gca,'FontWeight','bold','LineWidth',2,'FontSize',16)