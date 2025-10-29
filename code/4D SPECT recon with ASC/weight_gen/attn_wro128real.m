%attenuation weight generator
%Mingwu Jin, Aug. 01, 2006
%revised Sep. 23.
%No symmetry in attenuation case. Must 64 angles instead of 32 angles!

load weight128_120mn
load roi128
% dir_name='D:\imagereconstruction\4D\real_dataUM\P53606\';
% map_name='map_ostr.R35';
% fid=fopen([dir_name map_name],'rb');
% attn=fread(fid,128^3,'single');
% fclose(fid);
% attn=reshape(attn,[128 128 128]);
% save real_attnmap attn
load real_attnmap attn
x=attn(:,:,23:86);clear attn

    x=x.*repmat(roi,[1 1 64]);
    wp_attnwgt=cell(15360,1);%128*120
    tic
    for n=1:120
        if n<61
            wp_attnwgt((n-1)*128+1:n*128)=attnmat128(x,n,wp_vray((n-1)*128+1:n*128),...
                wp_ipxl((n-1)*128+1:n*128),wp_wgt((n-1)*128+1:n*128));
        else%No symmetry as in no attenuation case. Must have 64 angles instead of 32 angles!
            m=n-60;
            wp_attnwgt((n-1)*128+1:n*128)=attnmat128(x,n,wp_vray((m-1)*128+1:m*128),...
                wp_ipxl((m-1)*128+1:m*128),wp_wgt((m-1)*128+1:m*128));
        end;
    end,toc%15.7 seconds for 64 angles(local). %51 seconds for 120 angles (128*128*64)
    %normalize!!! Mar. 23,2006
    for j=1:120
        for i=1:128
            wp=wp_attnwgt{(j-1)*128+i};
            wp_S(i)=max(sum(wp));
        end
        wp_M(j)=max(wp_S);
    end%find peak for each angle
    for j=1:120
        for i=1:128
            wp_attnwgt{(j-1)*128+i}=wp_attnwgt{(j-1)*128+i}/max(wp_M);
        end
    end
    tic,%filename=['D:\imagereconstruction\4D\weights_gen\weight128_attn' num2str(g) '.mat'];
    filename='weight128_attn1.mat';
    save(filename,'wp_attnwgt');toc%38.7 seconds for saving 64 angles weights. %95sec for 120 angles (128*128*64)
%on the fly(more time/CPU demanding) or pre-stored(big files/memory demanding)!!