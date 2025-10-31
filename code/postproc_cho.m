function x = postproc_cho(im, s)

%post process saved 2D transversal slice of 4D images 
% crop and interpolates to s×s resolution for CHO analysis.

% im=im(3:end,:);
S = size(im); % current square! S(1)=S(2)

% Create meshgrid for interpolation
[x,y]=meshgrid(1:(S(2)-1)/(s-1):S(2),1:(S(1)-1)/(s-1):S(1));

% Interpolate to s×s
x=interp2(im,x,y);




