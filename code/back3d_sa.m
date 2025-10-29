function y=back3d_sa(pj,ang,wp_vray,wp_ipxl,wp_wgt,blur,gb_temp)
%3D backprojection with blurring at single angle
%64*64 --> 64*64*64
[S,M]=size(pj);%M: ray #; S: slice #
N=M;%plane
Total_stops=64;
%for second half backprojection
%may be not necessary if both proj3d_sa.m and back3d_sa.m do not include
%the mirror.
if ang>Total_stops/2 %revised Mingwu Jin, May 31, 2006
    pj=fliplr(pj);
end

ppj=zeros(S,M,N);
if blur==1
    for n=1:N
        ppj(:,:,n)=conv2(pj,gb_temp{n},'same');
    end
else
    ppj=repmat(pj,[1 1 N]);
end%blurred projection for each plane
ppj=permute(ppj,[2 1 3]);

%assign each plane to image coordinates
y=zeros(M*N,S);
for j=1:M%rays
    wp_v=wp_vray{j};
    if ~isempty(wp_v)
        if ang>Total_stops/2%for second half projection%revised Mingwu Jin, May 31, 2006
            wp_v=N-wp_v+1;
        end
        wp_i=wp_ipxl{j};
        wp_w=wp_wgt{j};
        L=length(wp_v);
        if size(wp_w,1)==1
            wp_w=repmat(wp_w',1,S);%without attenuation! %Mingwu Jin, Aug. 01, 2006
        end
        for i=1:L
            y(wp_i(i),:)=y(wp_i(i),:)+ppj(j,:,wp_v(i)).*wp_w(i,:);
        end %0.66 for notebook, 0.07 for cluster
    end
end
y=reshape(y,[M N S]);