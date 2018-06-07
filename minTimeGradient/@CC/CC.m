function  res = CC(csm)
% Operator to combine cols with a Roemer weighted reconstruction.
% Sum of squares can be performed with a only ones csm input.
% Handles 5D matrices with dimensions[x y z coils dynamics]
%
% Tom Bruijnen - University Medical Center Utrecht - 201609

res.S=csm;
res.adjoint=1; % 1 = forward (multicoil --> single coil). -1 = inverse operator
res=class(res,'CC');

%END
end