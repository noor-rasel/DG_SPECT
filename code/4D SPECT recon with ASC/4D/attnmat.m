function wp_attnwgt=attnmat(attn_map,ang,wp_vray,wp_ipxl,wp_wgt)
%attenuated transition matrix weights at angle "ang" (1-32)
%one slice weight --> 64 slices due to different attenuation maps.
%Mingwu Jin, Aug. 01, 2006

[M,N,S]=size(attn_map);%detector/ray # (x) = plane # (y) = slice # (z)
FOV=40.0849;%~=40.5=64*.634

attn_map=reshape(attn_map,M*N,S);
pj=zeros(M,S,M);%plane, slice, detector/ray
for j=1:M
    wp_v=wp_vray{j};   
    if ~isempty(wp_v)
        wp_i=wp_ipxl{j};
        wp_w=wp_wgt{j};
        L=length(wp_v);
        wp_w=repmat(wp_w',1,S);
        for i=1:L
            pj(wp_v(i),:,j)=pj(wp_v(i),:,j)+attn_map(wp_i(i),:).*wp_w(i,:);
        end %0.66 for notebook, 0.07 for cluster       
    end
end
%pj=permute(pj,[2 3 1]);%Order: S (slice), detector (M), plane (M)
%discretized integral
pj=cumsum(pj,1);
pj=exp(-pj*FOV);%normalized weights are scaled back
%get attenuated weights
wp_attnwgt=cell(64,1);
for j=1:M
    wp_v=wp_vray{j};   
    if ~isempty(wp_v)
        wp_i=wp_ipxl{j};
        wp_w=wp_wgt{j};
        L=length(wp_v);
        wp_w=repmat(wp_w',1,S);
        for i=1:L
            wp_w(i,:)=pj(wp_v(i),:,j).*wp_w(i,:);
        end %0.66 for notebook, 0.07 for cluster       
    end
    wp_attnwgt{j}=wp_w;%L*S matrix. L is varying.
end