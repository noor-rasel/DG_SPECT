function rim=fbp_3dMALL(nsin,px_num,s_ang,e_ang,cutoff)
%fbp 3D on MATLAB backprojection.

if nargin<5 cutoff=1;end
[M,N,f,s]=size(nsin);

tcount=sum(nsin(:));
%ramp filter
Np=max(N,2^nextpow2(N));
H = designFilter('ram-lak', Np, cutoff);
H=H./sum(H.^2);
H=repmat(H',M,1);

%filtering
for i=1:f
   for j=1:s      
      p = nsin(:,:,i,j);    % p holds fft of projections        
      if Np>N
          p(:,Np)=0;
      end% Zero pad projections 
%       p(M*2,Np*2)=0;
      p=fft(p,[],2);
      p=p.*H;
      p = real(ifft(p,[],2));     % p is the filtered projections
      nsin(:,:,i,j)=p(1:M,1:N);
   end
end

if px_num==0
    load weight64_mn
    pxn=64;
    tot_ang=64;
    rim=zeros(N,N,M,s);
elseif px_num==1
    load weight128_120mn
    pxn=128;
    tot_ang=120;
    rim=zeros(N,N,M,s);
end

nn=s_ang:e_ang;
for j=1:s
    fnsin=nsin(:,:,:,j);
    for sino_num=1:f; %f==(e_ang-s_ang+1)
        n=nn(sino_num);        
        if n<tot_ang/2+1
            rim(:,:,:,j)=rim(:,:,:,j)+back3d_saALL(fnsin(:,:,sino_num),n,wp_vray((n-1)*pxn+1:n*pxn),...
                wp_ipxl((n-1)*pxn+1:n*pxn),wp_wgt((n-1)*pxn+1:n*pxn),0,0,tot_ang);
        else
            m=n-tot_ang/2;
            rim(:,:,:,j)=rim(:,:,:,j)+back3d_saALL(fnsin(:,:,sino_num),n,wp_vray((m-1)*pxn+1:m*pxn),...
                wp_ipxl((m-1)*pxn+1:m*pxn),wp_wgt((m-1)*pxn+1:m*pxn),0,0,tot_ang);
        end
    end
end

% rim(rim<0)=0;
% rim=rim*tcount/sum(rim(:));