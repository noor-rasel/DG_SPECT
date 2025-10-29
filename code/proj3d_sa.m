function pj=proj3d_sa(x,ang,wp_vray,wp_ipxl,wp_wgt,blur,gb_temp)
%3D projection with blurring at single angle
%64*64*64 --> 64*64

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
%         wp_v=wp_v+1;
        wp_i=wp_ipxl{j};
        wp_w=wp_wgt{j};
        L=length(wp_v);
%         [ix,iy]=ind2sub([64 64],wp_i);
%         ix=ix+1;
%         wp_i=sub2ind([64 64],ix,iy);%slower
%         dc=64*ones(1,L);
%         iy=floor(wp_i/64);
%         ix=mod(wp_i,dc);ix=ix+1;
%         wp_i=iy*64+ix;%could adjust before loading!
        %method 1: matrix
%         temp=zeros(L*M,S);
%         ind=1:L;ind=ind+(wp_v-1)*L;%converting pixel coordinates to 1:L.
%         temp(ind,:)=x(wp_i,:).*repmat(wp_w',1,S);
%         temp=reshape(temp,[L M S]);
%         pj(:,:,j)=pj(:,:,j)+squeeze(sum(temp));%0.43 for cluster, much slower!
        %method 2: loop
        if size(wp_w,1)==1
            wp_w=repmat(wp_w',1,S);%No attenuation case! %Mingwu Jin, Aug. 01, 2006
        end
        for i=1:L
            pj(wp_v(i),:,j)=pj(wp_v(i),:,j)+x(wp_i(i),:).*wp_w(i,:);
        end %0.66 for notebook, 0.07 for cluster       
    end
end
pj=permute(pj,[2 3 1]);%Order: S (slice), detector (M), plane (M)
if blur==1
    for n=1:M
        pj(:,:,n)=conv2(pj(:,:,n),gb_temp{n},'same');%filter2(gb_temp{n},pj(:,:,n));
    end
end
pj=sum(pj,3);
%for second half projection
%may be not necessary if both proj3d_sa.m and back3d_sa.m do not include
%the mirror.
if ang>Total_stops/2%revised Mingwu Jin, May 31, 2006
    pj=fliplr(pj);
end