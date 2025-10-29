%attenuation weight generator
%Mingwu Jin, Aug. 01, 2006
%revised Sep. 23.
%No symmetry in attenuation case. Must 64 angles instead of 32 angles!

load weight64_mn
load roi
G=1;%16; 16 gates are similar. Approximation.
for g=1:G
    tic,filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'r');
    x=fread(fid,64^3,'single');
    fclose(fid);
    x=reshape(x,[64 64 64]);
    x=x.*repmat(roi,[1 1 64]);
    wp_attnwgt=cell(4096,1);%64*64
    for n=1:64
        if n<33
            wp_attnwgt((n-1)*64+1:n*64)=attnmat(x,n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64));
        else%No symmetry as in no attenuation case. Must have 64 angles instead of 32 angles!
            m=n-32;
            wp_attnwgt((n-1)*64+1:n*64)=attnmat(x,n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64));
        end;
    end,toc%15.7 seconds for 64 angles.
    %normalize!!! Mar. 23,2006
    for j=1:64
        for i=1:64
            wp=wp_attnwgt{(j-1)*64+i};
            wp_S(i)=max(sum(wp));
        end
        wp_M(j)=max(wp_S);
    end%find peak for each angle
    for j=1:64
        for i=1:64
            wp_attnwgt{(j-1)*64+i}=wp_attnwgt{(j-1)*64+i}/max(wp_M);
        end
    end
    tic,filename=['D:\imagereconstruction\4D\weights_gen\weight64_attn' num2str(g) '.mat'];
    save(filename,'wp_attnwgt');toc%38.7 seconds for saving 64 angles weights.
end
%on the fly(more time/CPU demanding) or pre-stored(big files/memory demanding)!!