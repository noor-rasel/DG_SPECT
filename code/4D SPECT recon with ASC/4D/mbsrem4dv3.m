function [x,v_num,obj_func]=mbsrem4dv3(sino,x,sub,iter,of_tag,sbeta,gbeta,blur,attn,scat,motion_tag,noisenum)
%Modified block sequential regularized EM.
%3D and 4D with motion compensation.
%Mingwu Jin, July 2006
%add attenuation correction (AC) and scatter correction (SC) Oct 07,2006.
%of_tag=1: used to calcualte MAP objective function and guarantee its
%increasing.
%single precision: Nov. 09, 2006
%add true motion: Nov. 13, 2006

if nargin<8
    blur=1;
    attn=0;
    scat=0;
    motion_tag=0;%0: estimated motion; 1: true NCAT motion.
end%default SPECT without AC and SC
x=single(x);
sino=single(sino);
scat=single(scat);
[nx,ny,nz,G]=size(x);%G gates.
[ns,nr,T,G]=size(sino);%nr=nx=ny(bin 64); ns=nz(slice 64); T (angle: 64).

%spatial prior template
if sbeta~=0
    %prior 1:
    sfilt=-ones(3,3,3);sfilt(2,2,2)=26;
    sfilt=sfilt/26;
    %prior 2:
    %     sfilt=zeros(3,3,3);
    %     sfilt(2,2,2)=6;
    %     sfilt([5 11 13 15 17 23])=-1;
end
%index for 5 gates Montion compensation
if gbeta~=0
    for g=1:G
        ind5(g,:)=mod((g-2:g+2),G);
    end
    ind5(ind5==0)=G;
    ind5(:,3)=[];
end


%generate subset index:"index";
%         number of the elements in the sub-set: "sub_com"
if sub==1
    index=1:T;index=index';
    sub_com=T;
else
    sub_com=T/sub;
    s_sub=reshape(1:(sub_com/4)*sub,sub_com/4,sub);
    index=zeros(sub_com/4,4,sub);
    for n=1:sub
        index(:,:,n)=repmat(s_sub(:,n),1,4)+repmat((0:3)*(sub_com/4)*sub,sub_com/4,1);
    end
    n1=1:sub/2;n2=sub/2+1:sub;
    n3=zeros(sub,1);
    for n=1:sub/2
        n3(2*n-1:2*n)=[n1(n) n2(n)];
    end
    index=index(:,:,n3);
    index=reshape(index,sub_com,sub);
    %sub=2,4,8,16 for 64 angles verified. Each column of "index" is the
    %index for one sub-set.
end

%modified BSREM global parameters
%(1)Upper bound of reconstructed voxel intensity:BSREM II.
load projMat_min
for g=1:G
    temp=sino(:,:,:,g)./pj;
    int_U(g)=max(temp(:));
end
int_U=max(int_U);clear pj temp
%(2)lower bound
epsilon=eps;
% nz_sino=find(sino>0);
% epsilon=min(sino(nz_sino));
% temp=ones(64,64,64);
% y=zeros(ns,nr,T,G);

%attenuation weights
% tic
if attn==1
    load weight64_mn
    load weight64_attn1.mat
    wp_wgt=wp_attnwgt;clear wp_attnwgt;
else
    load weight64_mn %no attenuation
end
% log_name=['wei_time_n' num2str(noisenum+1) '.mat'];
% wei_time=toc;save(log_name,'wei_time');
%load gbk64
% tic
load gbk64_umass
% log_name=['blu_time_n' num2str(noisenum+1) '.mat'];
% blu_time=toc;save(log_name,'blu_time');
% %projection of the identity 3D image
% for g=1:G
%     for n=1:T
%         if n<33
%             y(:,:,n,g)=proj3d_sa(x(:,:,:,g),n,wp_vray((n-1)*64+1:n*64),...
%                 wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
%         else
%             m=n-32;
%             y(:,:,n,g)=proj3d_sa(x(:,:,:,g),n,wp_vray((m-1)*64+1:m*64),...
%                 wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
%         end
%     end
% end
% nz_ind=find(y>0);
% %Spatial prior & temporal prior are both equal to 0 as the identity 3D
% %image used. No need for them to calculated objective function.
% LOF_iden=sum(sino(nz_ind).*log(y(nz_ind))-y(nz_ind));
% y=sino(nz_sino).*log(sino(nz_sino))-sino(nz_sino);
% LOF_iden=LOF_iden-sum(y(:));
% y=y+LOF_iden;
% y=exp(y(nz_sino)./sino(nz_sino));
% y=min(y);
% epsilon=min(y,epsilon);
% if epsilon==0%always satisfied
%     epsilon=eps;
% end
%(3)sub-seted denominator
rsc=zeros(nx,ny,nz,'single');
for n=1:64
    if n<33
        rsc=rsc+back3d_sa(ones(nz,nx,'single'),n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
    else
        m=n-32;
        if attn==1%attenuation weights
            rsc=rsc+back3d_sa(ones(nz,nx,'single'),n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
        else
            rsc=rsc+back3d_sa(ones(nz,nx,'single'),n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
        end
    end
end
zz_rsc=find(rsc==0);
rsc=rsc/sub;
%load motion compensation matrix: save time with the price of memory: 1GB
%extra memory needed. Save 16(g)*16(sub)*10(iter)*1.5(sec, load time)=3840 seconds

if gbeta~=0 & motion_tag==0
%     tic
    for g=1:G
        motionfile=['n4dMM' num2str(g) '_n' num2str(noisenum)];
        load(motionfile);
        M{g}=MM;
    end
%     log_name=['mot_time_n' num2str(noisenum+1) '.mat'];
%     mot_time=toc;save(log_name,'mot_time');
    clear MM
elseif gbeta~=0 & motion_tag==1
    for g=1:G
        motionfile=['true_4dMM' num2str(g)];
        load(motionfile);
        M{g}=MM;
    end
    clear MM
end

load roi
alpha0=1;
%objective function:time consuming. use for parameter determine.
%Afterwards not necessary! Necessary for adaptive alpha determine.
k=1;
if of_tag==1
    ps=zeros(G,1,'single');
    pg=zeros(G,1,'single');
    g_sub=zeros(ns,nr,T,G,'single');
    for g=1:G
        for m=1:T
            if m<33
                g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                    wp_vray((m-1)*64+1:m*64),...
                    wp_ipxl((m-1)*64+1:m*64),...
                    wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
            else
                comp_m=m-32;
                if attn==1%attenuation weights
                    g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                        wp_vray((comp_m-1)*64+1:comp_m*64),...
                        wp_ipxl((comp_m-1)*64+1:comp_m*64),...
                        wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
                else
                    g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                        wp_vray((comp_m-1)*64+1:comp_m*64),...
                        wp_ipxl((comp_m-1)*64+1:comp_m*64),...
                        wp_wgt((comp_m-1)*64+1:comp_m*64),blur,gb_temp);
                end
            end
        end
        if sbeta~=0
            %spatial prior
            %prior 1: (Xi-1/Ni*sum(Xj))^2, Xj is the neighbor f Xi
            temp=imfilter(x(:,:,:,g),sfilt,'replicate','same');
            temp=temp.*repmat(roi,[1 1 64]);
            ps(g)=sbeta*temp(:)'*temp(:);
            %prior 2: first oder quadratic. sum((Xi-Xj)^2)/6
            %                 neigh_index=[5 11 13 15 17 23];
            %                 for nn=1:6
            %                     temp_filt=zeros(3,3,3);
            %                     temp_filt(neigh_index(nn))=-1;
            %                     temp_filt(2,2,2)=1;
            %                     temp=imfilter(x(:,:,:,g),temp_filt,'replicate','same');
            %                     temp=temp.*repmat(roi,[1 1 64]);%dsp(temp)
            %                     ps(g)=ps(g)+sbeta*temp(:)'*temp(:)/6/2;
            %                 end
        else
            ps=0;
        end
        if gbeta~=0
            %temporal prior
            %             motionfile=['testMM' num2str(g)];
            %             load(motionfile);
            %             temp=x(:)'*MM;
            %             pg(g)=gbeta*temp*temp';
            temp=double(x(:,:,:,ind5(g,:)));
            temp=M{g}*temp(:);
            temp1=x(:,:,:,g);temp=temp1(:)-temp(:);
            pg(g)=gbeta*temp'*temp;
        else
            pg=0;
        end
    end
    %scatter compensation
    g_sub=g_sub+scat;
    nz_g=find(g_sub>0);
    obj_func(k)=sum(sino(nz_g).*log(g_sub(nz_g))-g_sub(nz_g))-sum(ps)-sum(pg);
else
    obj_func=[];
end
%main loop
v_num=0;%monitor violation numbers: n=v_num+k;
for n=1:iter
    alpha=alpha0/(k+1);%0.015/n^.33*2;
    %most important for convergence! Need to adjust for different case.
    %1/(n/15+1) for 8 subsets in Ahn paper, good for EM 16 subsets 10 iter.
    %1/(n*4+1) for sbeta=0.01. 5 iterations
    if of_tag==1
        x_buff=x;
    end
    for s=1:sub
        y=zeros(nx,ny,nz,G,'single');
        for g=1:G
            temp=zeros(nx,ny,nz,'single');
            for m=1:sub_com
                %sino_sub=zeros(nz,nx);
                g_sub=zeros(nz,nx,'single');
                sino_sub=sino(:,:,index(m,s),g);
                if index(m,s)<33
                    g_sub=proj3d_sa(x(:,:,:,g),index(m,s),...
                        wp_vray((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_ipxl((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                else
                    comp_ind=index(m,s)-32;
                    if attn==1%attenuation weights
                        g_sub=proj3d_sa(x(:,:,:,g),index(m,s),...
                            wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                            wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                            wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                    else
                        g_sub=proj3d_sa(x(:,:,:,g),index(m,s),...
                            wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                            wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                            wp_wgt((comp_ind-1)*64+1:comp_ind*64),blur,gb_temp);
                    end
                end
                %scatter compensation
                if length(scat)>1
                    g_sub=g_sub+scat(:,:,index(m,s),g);
                end
                %g_sub(i)>epsilon
                nnz=find(g_sub>epsilon);
                %sino(i)>0 and g_sub(i)<=epsilon
                nnz_small=find(g_sub<=epsilon&sino_sub~=0);%seems never to happen
                g_sub(nnz)=sino_sub(nnz)./g_sub(nnz)-1;
                if ~isempty(nnz_small)
                    g_sub(nnz_small)=(sino_sub(nnz_small)/epsilon-1)...
                        -sino_sub(nnz_small)/epsilon^2.*(g_sub(nnz_small)-epsilon);
                end
                %                 z_sino=find(sino_sub==0);
                %                 if length(nnz)<(nz*nx-length(z_sino))%seems never to happen
                %                     nnz=find(g_sub<=epsilon);
                %                     nnz(z_sino)=[];
                %                     g_sub(nnz)=(sino_sub(nnz)/epsilon-1)...
                %                         -sino_sub(nnz)/epsilon^2*(g_sub(nnz)-epsilon);
                %                 end
                %derivative of objective function:backprojection of ratio
                if index(m,s)<33
                    temp=temp+back3d_sa(g_sub,index(m,s),...
                        wp_vray((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_ipxl((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                else
                    comp_ind=index(m,s)-32;
                    if attn==1%attenuation weights
                        temp=temp+back3d_sa(g_sub,index(m,s),...
                            wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                            wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                            wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                    else
                        temp=temp+back3d_sa(g_sub,index(m,s),...
                            wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                            wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                            wp_wgt((comp_ind-1)*64+1:comp_ind*64),blur,gb_temp);
                    end
                end
            end
            %derivative of objective function:backprojection of ratio
            y(:,:,:,g)=temp;
            if sbeta~=0
                %derivative of objective function:spatial prior
                %prior 1:
                temp=imfilter(x(:,:,:,g),sfilt,'replicate','same');
                temp=2*sbeta*imfilter(temp,sfilt,'replicate','same');
                %prior 2:
                %temp=2*sbeta*imfilter(x(:,:,:,g),sfilt,'replicate','same')/6;
                y(:,:,:,g)=y(:,:,:,g)-temp/sub;
            end
            if gbeta~=0 
                %derivative of objective function:temporal prior (5 gates)
%                 motionfile=['testMM' num2str(g)];
%                 load(motionfile);
%                 temp=2*gbeta*x(:)'*MM;temp=reshape(temp,[64,64,64]);
%                 y(:,:,:,g)=y(:,:,:,g)-temp/sub;
                temp=double(x(:,:,:,ind5(g,:)));
                temp=M{g}*temp(:);temp1=x(:,:,:,g);
                temp=2*gbeta*(temp1(:)-temp(:));
                temp=reshape(temp,[nx ny nz]);                
                y(:,:,:,g)=y(:,:,:,g)-temp/sub;
            end
            %step size of derivative
            temp=x(:,:,:,g);
            temp1=zeros(nx,ny,nz,'single');
            upper_part=find(temp>=int_U/2);
            temp1(upper_part)=(int_U-temp(upper_part))./rsc(upper_part);
            lower_part=find(temp<int_U/2&rsc>0);
            temp1(lower_part)=temp(lower_part)./rsc(lower_part);
            y(:,:,:,g)=y(:,:,:,g).*temp1;
        end
        %update
        x=x+alpha*y;
        x(x<=0)=eps;%may apply ROI
        x(x>=int_U)=int_U-eps;
    end
    if of_tag==1%to make sure increasing of objective function
        ps=zeros(G,1,'single');
        pg=zeros(G,1,'single');
        g_sub=zeros(ns,nr,T,G,'single');
        for g=1:G
            for m=1:T
                if m<33
                    g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                        wp_vray((m-1)*64+1:m*64),...
                        wp_ipxl((m-1)*64+1:m*64),...
                        wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
                else
                    comp_m=m-32;
                    if attn==1%attenuation weights
                        g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                            wp_vray((comp_m-1)*64+1:comp_m*64),...
                            wp_ipxl((comp_m-1)*64+1:comp_m*64),...
                            wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
                    else
                        g_sub(:,:,m,g)=proj3d_sa(x(:,:,:,g),m,...
                            wp_vray((comp_m-1)*64+1:comp_m*64),...
                            wp_ipxl((comp_m-1)*64+1:comp_m*64),...
                            wp_wgt((comp_m-1)*64+1:comp_m*64),blur,gb_temp);
                    end
                end
            end
            if sbeta~=0
                %spatial prior
                %prior 1: (Xi-1/Ni*sum(Xj))^2, Xj is the neighbor f Xi
                temp=imfilter(x(:,:,:,g),sfilt,'replicate','same');
                temp=temp.*repmat(roi,[1 1 64]);
                ps(g)=sbeta*temp(:)'*temp(:);
                %prior 2: first oder quadratic. sum((Xi-Xj)^2)/6
                %             neigh_index=[5 11 13 15 17 23];
                %             for nn=1:6
                %                 temp_filt=zeros(3,3,3);
                %                 temp_filt(neigh_index(nn))=-1;
                %                 temp_filt(2,2,2)=1;
                %                 temp=imfilter(x(:,:,:,g),temp_filt,'replicate','same');
                %                 temp=temp.*repmat(roi,[1 1 64]);%dsp(temp)
                %                 ps(g)=ps(g)+sbeta*temp(:)'*temp(:)/6/2;
                %             end
            else
                ps=0;
            end
            if gbeta~=0
                %temporal prior
                %             motionfile=['testMM' num2str(g)];
                %             load(motionfile);
                %             temp=x(:)'*MM;
                %             pg(g)=gbeta*temp*temp';
                temp=double(x(:,:,:,ind5(g,:)));
                temp=M{g}*temp(:);
                temp1=x(:,:,:,g);temp=temp1(:)-temp(:);
                pg(g)=gbeta*temp'*temp;
            else
                pg=0;
            end
        end
        %scatter compensation
        g_sub=g_sub+scat;
        nz_g=find(g_sub>0);
        obj_buff=sum(sino(nz_g).*log(g_sub(nz_g))-g_sub(nz_g))-sum(ps)-sum(pg);
        if obj_buff<obj_func(k)
            alpha0=alpha0/2;
            x=x_buff;
            v_num=v_num+1;
        else
            k=k+1;
            obj_func(k)=obj_buff;
        end
    else
        k=k+1;
    end
end