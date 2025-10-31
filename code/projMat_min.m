function pj=projMat_min(ang,wp_vray,wp_ipxl,wp_wgt,blur,gb_temp)
%Find minimal(>0) of 3D projection matrix at "ang"(1-64).
%Used for determining the upper bound for BSREMII.
%blur==1, SPECT(distance-dependent blur); else, PET.
%
%Mingwu Jin, July 2006.

x=ones(64,64,64);
[M,N,S]=size(x);%detector/ray # (x) = plane # (y) = slice # (z)
Total_stops=64;%revised Mingwu Jin, May 31, 2006

x=reshape(x,M*N,S);
pj=zeros(M,S,M);%plane, slice, detector/ray
for j=1:M
    wp_v=wp_vray{j};   
    if ~isempty(wp_v)
        if ang>Total_stops/2%for second half projection:%revised Mingwu Jin, May 31, 2006
            wp_v=N-wp_v+1;
        end
        wp_i=wp_ipxl{j};
        wp_w=wp_wgt{j};
        L=length(wp_v);
        wp_w=repmat(wp_w',1,S);
        for i=1:L
            pj(wp_v(i),:,j)=pj(wp_v(i),:,j)+x(wp_i(i),:).*wp_w(i,:);
        end %0.66 for notebook, 0.07 for cluster       
    end
end
pj=permute(pj,[2 3 1]);%Order: S (slice), detector (M), plane (M)
if blur==1
    for n=1:M
        pj(:,:,n)=conv2(pj(:,:,n),gb_temp{n},'same');
    end
end
pj(pj==0)=max(pj(:));
pj=min(pj,[],3);
%for second half projection
%may be not necessary if both proj3d_sa.m and back3d_sa.m do not include
%the mirror.
if ang>Total_stops/2%revised Mingwu Jin, May 31, 2006
    pj=fliplr(pj);
end