function rim=fbp_3dM(nsin,cutoff)
%fbp 3D on MATLAB backprojection.

if nargin<2 cutoff=1;end
[M,N,f,s]=size(nsin);
rim=zeros(M,N,f,s);
tcount=sum(nsin(:));
%ramp filter
Np=max(N,2^nextpow2(N));
H = designFilter('ram-lak', Np, cutoff);
H=H./sum(H.^2);
% H=repmat(H',M,1);
H=repmat(H',N,1);

%filtering
for i=1:f
   for j=1:s      
      p = nsin(:,:,i,j);    % p holds fft of projections        
      if Np>N
          p(:,Np)=0;
      end% Zero pad projections 
%       p(M*2,Np*2)=0;
      p=fft(p,[],2);
      p=p.*H';
%       p=p.*H;
      p = real(ifft(p,[],2));     % p is the filtered projections
      nsin(:,:,i,j)=p(1:M,1:N);
   end
end

load weight64_mn
for j=1:s
    fnsin=nsin(:,:,:,j);
        for n=1:64
        if n<33
            rim(:,:,:,j)=rim(:,:,:,j)+back3d_sa(fnsin(:,:,n),n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),0,0);
        else
            m=n-32;
            rim(:,:,:,j)=rim(:,:,:,j)+back3d_sa(fnsin(:,:,n),n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),0,0);
        end
    end
end

% rim(rim<0)=0;
% rim=rim*tcount/sum(rim(:));