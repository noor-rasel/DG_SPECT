% Erwan Gravier
% Research project
% 11/26/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the 3D butterworth
% filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage:  f=butter3d(sze,cutoff2,n)
% 
% where: sze    is the size of the image (sze*sze*sze)
%        cutoff is the cutoff frequency of the filter 0 - 1
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%        f: the coefficient of the filter to be multiplied by the image.
%        the frequency center is the center of the cube... (64,64,64)
% Tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=butter3d(sze,cutoff,n)
  
if cutoff<0 | cutoff>1
    error('cutoff frequency must be between 0 and 1');
end
    
% if rem(n,1) ~= 0 | n < 1
%     error('n must be an integer >= 1');
% end% Mingwu Jin, July 06, 2006

step=2/(sze);
% X, Y and Z matrices with ranges normalised to +/- 1
[X,Y,Z]=meshgrid(-1:step:1-step,-1:step:1-step,-1:step:1-step);

% matrix with the value of the radius for each pixel (origin at the center
% of the cube)
radius=sqrt(X.^2 +Y.^2+Z.^2);

% The filter
f=1./(1+(radius./cutoff).^(2*n));

