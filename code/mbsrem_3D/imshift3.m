function imr = imshift3(im,vect)
%% shift image according the vectors in 3D
%% linear interpolation.




maxshift = round(max(abs(vect)))+2;

[nx,ny,nz,G] = size(im);



imtemp = zeros(nx+2*maxshift,ny+2*maxshift,nz+2*maxshift);

imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift) = im;

imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift) = ...
    abs(vect(1)-(floor(vect(1))+1))*imtemp(floor(1+vect(1):nx+vect(1))+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift)...
    +abs(vect(1)-floor(vect(1)))*imtemp(floor(1+vect(1):nx+vect(1))+1+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift);

imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift) = ...
    abs(vect(2)-(floor(vect(2))+1))*imtemp(1+maxshift:nx+maxshift,floor(1+vect(2):ny+vect(2))+maxshift,1+maxshift:nz+maxshift)...
    +abs(vect(2)-floor(vect(2)))*imtemp(1+maxshift:nx+maxshift,floor(1+vect(2):ny+vect(2))+1+maxshift,1+maxshift:nz+maxshift);

imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift) = ...
    abs(vect(3)-(floor(vect(3))+1))*imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,floor(1+vect(3):nz+vect(3))+maxshift)...
    +abs(vect(3)-floor(vect(3)))*imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,floor(1+vect(3):nz+vect(3))+1+maxshift);

imr = imtemp(1+maxshift:nx+maxshift,1+maxshift:ny+maxshift,1+maxshift:nz+maxshift);

