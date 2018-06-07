function [I,Ist] = generate_cartesian_acquisition(F,W,B0,der_B0,ph)
% See what happens in the cartesian acquisition

% Rotate fat and non-fat in and out of phase
fat_phi=@(tau,woff)(sin((ph.woff_fat+woff).*tau));
water_phi=@(tau,woff)(sin(woff.*tau));

% Standard Cartesian acquisition
k=zeros(ph.N,ph.N);
for fe=1:ph.N
    t=ph.te+(fe-1)*ph.adc_dt;
    % Apply phase ramping due to off res + chem shift
    x=exp(-1j*2*pi*fat_phi(t,B0)).*F+...
        exp(-1j*2*pi*water_phi(t,B0)).*W;
    
    % Apply intravoxel dephasing 
    for xcoord=1:ph.N^2
    spins=t*der_B0(xcoord);
    if spins(end)>1
        x(xcoord)=0;
    else
        %x(xcoord)=x(xcoord)*exp(-1.6*abs(spins));
    end
    end
    
    % Store image per time-point
    Ist(:,:,1,1,fe)=x;
    
    % Sample a complete phase-encode line (say TE)
    X=ifftshift(fft2(fftshift(x)));
    k(fe,:)=X(fe,:);
end

% Make it more realistic by introducing gibbs ringing
k([1 end],:)=0;k(:,[1 end])=0;

% Fourier transform
I=ifftshift(ifft2(fftshift(k)));

if ph.vis==1
    subplot(4,4,[11 12 15 16]);imshow(abs(I),[]);title('Cartesian recon');
end

% END
end