function [W, F] = multi_compartmentalise(img,varargin)
% f is fat fraction

if nargin > 1
    f=varargin{1};
else
    f=.2;
end

W=img;
F=img;
W(W>36 & W<37)=f*36.6423187;
F(F<36 | F>37)=0;
F(F>36 & F<37)=(1-f)*36.6423187;

% END
end