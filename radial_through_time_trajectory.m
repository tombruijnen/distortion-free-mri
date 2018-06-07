function [k,dcf] = radial_through_time_trajectory(ph)
% Extract radial k-space coordinates, this does not include the generation
% of gradient waveforms.

% Preallocate k [ns nl 1 1 nt]
k=zeros(ph.N,round(ph.N/ph.R),1,1,ph.N);

% Golden increment between time-points
ga=pi/((1+sqrt(5))/2);

% Uniform increment within time-point
u=2*pi/round(ph.N/ph.R);

% Now compute all the sample points on single spoke
kx=linspace(0,ph.N-1,ph.N)'-(ph.N-1)/2;

% Modulate the phase of all the successive spokes
for dyn=1:ph.N
    for l=1:round(ph.N/ph.R)
        ang=(dyn-1)*ga+(l-1)*u;
        k(:,l,1,1,dyn)=kx*exp(-1j*ang);
    end
end

% Normalise to reconstruction matrix
k=ph.N*k/max(abs(k(:)));

% Density function
dcf=abs(k);dcf(floor(ph.N/2)+1,:)=1/ph.N;

% Split in [3 ns nl 1 1 ndyn]
knew=zeros([3,size(k)]);
knew(1,:,:,:,:,:)=real(k);knew(2,:,:,:,:,:)=imag(k);
k=knew;

% END
end