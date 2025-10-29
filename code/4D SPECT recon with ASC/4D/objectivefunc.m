function obj_func=objectivefunc(x,sino,blur,sbeta,gbeta)

[nx,ny,nz,G]=size(x);%G gates.
[ns,nr,T,G]=size(sino);
%spatial prior template
if sbeta~=0
    sfilt=-ones(3,3,3);sfilt(2,2,2)=26;
    sfilt=sfilt/26;
end
%index for 5 gates Montion compensation
if gbeta~=0
    for g=1:G
        ind5(g,:)=mod((g-2:g+2),G);
    end
    ind5(ind5==0)=G;
end

load weight64_mn
load gbk64

for g=1:G
    ps=zeros(G,1);
    pg=zeros(G,1);
    for m=1:T
        if m<33
            g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),...
                wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
        else
            comp_m=m-32;
            g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                wp_vray((comp_m-1)*64+1:comp_m*64),...
                wp_ipxl((comp_m-1)*64+1:comp_m*64),...
                wp_wgt((comp_m-1)*64+1:comp_m*64),blur,gb_temp);
        end
    end
    if sbeta~=0
        %spatial prior
%         for neigh_index=1:27
%             temp_filt=zeros(3,3,3);
%             if neigh_index~=14
%                 temp_filt(neigh_index)=sfilt(neigh_index);
%                 temp_filt(2,2,2)=sfilt(2,2,2)/26;
%                 temp=convn(x(:,:,:,g),temp_filt,'same');
%                 ps(g)=ps(g)+sbeta*temp(:)'*temp(:);
%             end
%         end%need to multiply 26
        temp=convn(x(:,:,:,g),sfilt,'same');
        temp=temp.*(repmat(roi,[1,1,64]));
        temp(:,:,1)=0;
        temp(:,:,64)=0;
        temp=temp.*x(:,:,:,g);
        ps(g)=sbeta*sum(temp(:));
    else
        ps=0;
    end
    if gbeta~=0
        %temporal prior
        temp=x(:,:,:,ind5(g,:));
        motionfile=['MM4dgate' num2str(g)];
        load(motionfile);
        temp=MM*temp(:);
        pg(g)=gbeta*temp'*temp;
    else
        pg=0;
    end
end
nz_g=find(g_sub>0);
obj_func=sum(sino(nz_g).*log(g_sub(nz_g))-g_sub(nz_g))-sum(ps)-sum(pg);