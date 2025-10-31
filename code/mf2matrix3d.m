function MM=mf2matrix3d(globalmot,g,gamma)

[M,N,S,D,G]=size(globalmot);         % size of the motion field
%triangular filter weights
k=1:G;
wei=abs(1-2*abs(g-k)/G).^gamma;
wei(g)=0;wei=wei/sum(wei);
%find circular shift index
ind=(1:G)';
ind=circshift(ind,G/2-g);
%left part: 1~G/2-1
tmov=zeros(M,N,S,3);
MM=cell(8);
M1=spalloc(M*N*S,M*N*S,8*M*N*S);
for j=1:G/2-1
    tmov=sum(globalmot(:,:,:,:,ind(j:G/2-1)),5);
    M1=motionmatrix3d(tmov(:,:,:,1),tmov(:,:,:,2),tmov(:,:,:,3));
    M1=wei(ind(j))*M1;
    MM{ind(j)}=M1';
%     MM(:,M*N*S*(ind(j)-1)+1:M*N*S*ind(j))=M1; %very slow: 451.008957seconds   
end
%right part: G/2~G-1
tmov=zeros(M,N,S,3);
for j=G/2:G-1
    tmov=sum(globalmot(:,:,:,:,ind(G/2:j)),5);
    M1=motionmatrix3d(-tmov(:,:,:,1),-tmov(:,:,:,2),-tmov(:,:,:,3));
    M1=wei(ind(j+1))*M1; 
    MM{ind(j+1)}=M1';
%     MM(:,M*N*S*(ind(j+1)-1)+1:M*N*S*ind(j+1))=wei(ind(j+1))*M1;    
end
MM{g}=spalloc(M*N*S,M*N*S,0);
%modified to generate 'I-M'
%MM{g}=-speye(M*N*S);
MM=cat(2,MM{:});
%MM=-MM;