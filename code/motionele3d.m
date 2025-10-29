% Erwan Gravier
% Research project 11/23/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is the implementation of the 3d optical flow
% algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vx_allT,vy_allT,vz_allT]=motionele3d(Sim,numb,lambda)
%                   Sim(:,:,:,1)--> reference frame (64*64*64)
%                   Sim(:,:,:,2)--> target frame    (64*64*64)
%                   numb: number of iterations
%                   lambda: smoothness parameter
%                   vi_allT: motion in the i direction (64*64*64)
%        ----->x
%        ..
%        .  .
%        .   .
%        y    z
% Time required to run the algorithm on 2 64*64*64 images:
% 51 iterations = 63.2s
% Tested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vx_allT,vy_allT,vz_allT]=motionele3d(Sim,numb,lambda)

% constants and variables
[M,N,L,K]=size(Sim);
vx_allT=zeros(M,N,L);     % normally 64*64*64
vy_allT=zeros(M,N,L);
vz_allT=zeros(M,N,L);
%normalization
Sim=Sim/(max(Sim(:)));

% kernel filters
Lx=zeros(3,3,3);
Ly=zeros(3,3,3);
Lz=zeros(3,3,3);
Lx(:,:,2)=[0 0 0;0 -1 1 ;0 -1 1];
Lx(:,:,3)=[0 0 0;0 -1 1 ;0 0 0];
Ly(:,:,2)=[0 0 0; 0 -1 -1; 0 1 1];
Ly(:,:,3)=[0 0 0; 0 -1 0; 0 1 0];
Lz(:,:,2)=[0 0 0; 0 -1 -1; 0 -1 0];
Lz(:,:,3)=[0 0 0; 0 1 1; 0 1 0];
% IMPORTANT-->because of the convolution we flip the kernels
Lx=flipdim(Lx,1);
Lx=flipdim(Lx,2);
Lx=flipdim(Lx,3);
Ly=flipdim(Ly,1);
Ly=flipdim(Ly,2);
Ly=flipdim(Ly,3);
Lz=flipdim(Lz,1);
Lz=flipdim(Lz,2);
Lz=flipdim(Lz,3);

% that one is symetrical, so no flipping needed...
Avg=zeros(3,3,3);
Avg(:,:,1)=[0 1/24 0; 1/24 1/12 1/24; 0 1/24 0];
Avg(:,:,2)=[1/24 1/12 1/24; 1/12 0 1/12; 1/24 1/12 1/24];
Avg(:,:,3)=[0 1/24 0; 1/24 1/12 1/24; 0 1/24 0];

tder=zeros(3,3,3);
tder(:,:,2)=[0 0 0; 0 1 1; 0 1 1];
tder(:,:,3)=[0 0 0; 0 1 1; 0 1 0];
tder=flipdim(tder,1);
tder=flipdim(tder,2);
tder=flipdim(tder,3);

% target, reference frames
tf=Sim(:,:,:,2);    % target frame (time k+1).
rf=Sim(:,:,:,1);    % reference frame (time k).

% spatial and temporal derivatives
Dx=1/6*(convn(rf,Lx,'same')+convn(tf,Lx,'same'));       % spatial derivative along x.
Dy=1/6*(convn(rf,Ly,'same')+convn(tf,Ly,'same'));       % spatial derivative along y.
Dz=1/6*(convn(rf,Lz,'same')+convn(tf,Lz,'same'));       % spatial derivative along y.
Dt=1/7*(convn(tf,tder,'same')-convn(rf,tder,'same'));    %temporal der.

% fprintf('preprocessing finished\n');

h=1;     %clock
while(1)
    vh1=convn(vx_allT,Avg,'same');  % calculate the weighted averages...
    vh2=convn(vy_allT,Avg,'same');
    vh3=convn(vz_allT,Avg,'same');
    % common and individual updates
    cupd=(Dx.*vh1+Dy.*vh2+Dz.*vh3+Dt)./(lambda.^2+Dx.^2+Dy.^2+Dz.^2);
    vx_allT=vh1-Dx.*cupd;
    vy_allT=vh2-Dy.*cupd;
    vz_allT=vh3-Dz.*cupd;

    % if enough iterations break...
    if h>=numb
    break;
    end
      
    h=h+1;
%     fprintf('iteration# %d done \n',h);
end
