%% Lets define and validate the operators

% Sensitivity operator
A=rand(200,200,1,8,1,5);
S=CC(randn(200,200,1,8));
SA=S*A;                    
A2=S'*A;                    % Works

% Fourier operator
I=repmat(phantom(128),[1 1 1 7 5]);FI=fftshift(fft2(fftshift(I)));
kx=linspace(-pi,pi,129);kx(end)=[];ky=linspace(-pi,pi,129);ky(end)=[];
[KX,KY]=meshgrid(kx,ky);
k=[];k(1,:,:,:,:,:)=KX;k(2,:,:,:,:,:)=KY;k(3,:,:,:,:,:)=zeros(size(KX));
k=repmat(k,[1 1 1 1 5]);
Id=[200 200 1 8 5 1 1 1 1 1 1 1];Kd=[128 128 1 1 5 1 1 1 1 1 1 1];
F=GG2D(k,Id,Kd);

A=F'*I;