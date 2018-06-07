function csm = generate_b0(nc,Rdims,varargin)
% Function to generate Gaussian profile based coil sensitivity maps
% uniformly distributed around a 2D image for an arbitrary number of coils.
% Only works for square images.
%
% 20160903 - University Medical Center Utrecht - Tom Bruijnen

% Default scale (sig)
if numel(varargin)<1
    sig=Rdims(1)/2;
else
    sig=varargin{1}*Rdims/10;
end

% Create gaussian formula
halfsize=Rdims(1)/2;
gaussian=@(x,y,dx,dy)(1/(sqrt(2*pi)*sig)*exp(-1*((x-dx-halfsize)^2+(y-dy-halfsize)^2)/(2*sig^2)));

% Circular uniform distributed
if numel(varargin)<2 
    % Calculate position of element based on angular distance
    dang=2*pi/nc;
    angles=0:dang:2*pi-dang;

    % Found corresponding position in img
    for c=1:nc
        dy1(c)=Rdims(1)/2*cos(angles(c));
        dx1(c)=Rdims(1)/2*sin(angles(c));

        dx(c)=(Rdims(1)/2)/sqrt(dx1(c)^2+dy1(c)^2)*dx1(c);
        dy(c)=(Rdims(1)/2)/sqrt(dx1(c)^2+dy1(c)^2)*dy1(c);
    end    
% Anterior coil abdominal style - only possible for n^2 coils
elseif strcmpi(varargin{2},'anterior') && mod(sqrt(nc),1)==0 
        dy=repmat(-Rdims(1)/2:Rdims(1)/(sqrt(nc)-1):Rdims(1)/2,[1 sqrt(nc)]);
        dx=repmat(-Rdims(1)/2:Rdims(1)/(sqrt(nc)-1):Rdims(1)/2,[sqrt(nc) 1]);
        dx=dx(:);dy=dy(:);
else
    fprintf('Error: No legit input ---> no output generated\n')
    csm=0;
    return;
end

% Generate sensitivity maps (translated Gaussians)
csm=zeros([Rdims nc]);
for c=1:nc
    for x=1:Rdims(1)    
        for y=1:Rdims(1)
            csm(x,y,c)=gaussian(x,y,dx(c),dy(c));
        end
    end
end

% Normalize and permute to [x y 1 coils]
csm=permute(csm/max(csm(:)),[1 2 3 4]);
%gtot=sum(csm,4);imshow(gtot,[])

end