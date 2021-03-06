function [F,W,B0,der_B0] = create_susceptibillity_phantom(W,F,ph,B0)
% Create susceptibillity phantom across sampling time domain

N=size(W,1); 

% Add susceptibillity
if ph.width_sus > 0
    dChi=zeros(N,N,N);
    dChi(N/2-ph.width_sus:N/2+ph.width_sus,N/2-ph.width_sus:N/2+...
        ph.width_sus,N/2-ph.width_sus:N/2+ph.width_sus)=3*10^-6+...
        .25*10^-6*rand(2*ph.width_sus+1,2*ph.width_sus+1,2*ph.width_sus+1); % Titanium alloy
    dB=ph.B0*calculateFieldShift(dChi,ph.voxel_size);
    dOhmega=dB*ph.gamma;
else
    dOhmega=zeros(N,N);
end

% If no B0 map is provided add linear B0 ramp
if nargin < 4
    B0=rot90(ph.off_res_max*normalise(generate_b0(1,N)),1);
end
    
% Sum them together
B0=B0+dOhmega(:,:,N/2);

% Calculate derivative of B0 for intravoxel dephasing
der_B0=1/4*conv2(B0,[0 -1 0; -1 4 -1; 0 -1 0],'same');

% END
end