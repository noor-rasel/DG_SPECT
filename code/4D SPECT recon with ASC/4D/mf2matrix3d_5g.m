function MM=mf2matrix3d_5g(globalmot,gamma)
%new version from mf2matrix3d.m(old one supports 8-gate 64*64*32, even number gates).
%Use 5-gate motion, smaller MM, g=3 (centered)(more general: oder number
%gates).
%Mingwu Jin, Oct 08, 2006

g=3;%the frame to be compensated is centered!!!
[M,N,S,D,G]=size(globalmot);         % size of the motion field
%triangular filter weights
k=1:G;
wei=abs(1-2*abs(g-k)/G).^gamma;
wei(g)=[];wei=wei/sum(wei);
%find circular shift index: to center the frame to be compensated. not
%neccessary.
% ind=(1:G)';
% ind=circshift(ind,G/2-g);
%15 16 1 2 3=>15-16 16-1 1-2 2-3 3-4: only first 4 motion needed! "frame 4"
%is not used.
%left part: 1~(G+1)/2-1
tmov=zeros(M,N,S,3);
MM=cell(4);
%M1=spalloc(M*N*S,M*N*S,8*M*N*S);
for j=1:(G+1)/2-1
    tmov=sum(globalmot(:,:,:,:,j:(G+1)/2-1),5);
    M1=motionmatrix3d(tmov(:,:,:,1),tmov(:,:,:,2),tmov(:,:,:,3));
    M1=wei(j)*M1;
    MM{j}=M1';
end
%right part: (G+1)/2~G-1
tmov=zeros(M,N,S,3);
for j=(G+1)/2:G-1
    tmov=sum(globalmot(:,:,:,:,(G+1)/2:j),5);
    M1=motionmatrix3d(-tmov(:,:,:,1),-tmov(:,:,:,2),-tmov(:,:,:,3));
    M1=wei(j)*M1; 
    MM{j}=M1';
end
%MM{g}=spalloc(M*N*S,M*N*S,0);
%modified to generate 'I-M'
%MM{g}=-speye(M*N*S);
MM=cat(2,MM{:});
%MM=-MM;