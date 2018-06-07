function [F,W,B0,der_B0] = create_susceptibillity_phantom(ph)
% Create susceptibillity phantom across sampling time domain

% Create water+fat layer
W=phantom(ph.N);F=zeros(size(W));
F(W==1)=ph.ff;W(W==1)=ph.wf; 
 
% Add susceptibillity
dChi=zeros(ph.N,ph.N,ph.N);
dChi(ph.N/2-ph.width_sus:ph.N/2+ph.width_sus,ph.N/2-ph.width_sus:ph.N/2+...
    ph.width_sus,ph.N/2-ph.width_sus:ph.N/2+ph.width_sus)=3*10^-6+...
    .25*10^-6*rand(2*ph.width_sus+1,2*ph.width_sus+1,2*ph.width_sus+1); % Titanium alloy
dB=ph.B0*calculateFieldShift(dChi,ph.voxel_size);
dOhmega=dB*ph.gamma;

% Add off-resonance field
B0=rot90(ph.off_res_max*normalise(generate_b0(1,ph.N)),1)+dOhmega(:,:,ph.N/2);

% Calculate derivative of B0 for intravoxel dephasing
der_B0=1/4*conv2(B0,[0 -1 0; -1 4 -1; 0 -1 0],'same');

% Visualization
if ph.vis==1
    subplot(4,4,[1 2 5 6]);imshow(W,[]);title('Water phantom')
    subplot(4,4,[3 4 7 8]);imshow(F,[]);title('Fat phantom')
    subplot(4,4,[9 10 13 14]);imshow(B0,[]);colormap(gca, jet);title('B0 - susceptibillity');
end
% END
end