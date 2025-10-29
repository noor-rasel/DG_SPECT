function [x,v_num,obj_func]=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat)
%Modified block sequential regularized EM.
%3D and 4D with motion compensation.
%Mingwu Jin, July 2006
%add attenuation correction (AC) and scatter correction (SC) Oct 07,2006.
%of_tag=1: used to calcualte MAP objective function and guarantee its
%increasing.
%single precision: Nov. 09, 2006
%add true motion: Nov. 13, 2006
%dedicated for IRIX clincal data

sub=17;%total 68 views divided into 17 sub-sets.

if nargin<7
    blur=1;
    attn=0;
    scat=0;
end%default SPECT without AC and SC
x=single(x);
sino=single(sino);
scat=single(scat);
[nx,ny,nz,G]=size(x);%G gates.
[ns,nr,T,G]=size(sino);%nr=nx=ny(bin 128); ns=nz(slice 64); T (angle: 68).

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
    s_sub=[1 10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 9];
    index=[s_sub;s_sub+sub;s_sub+sub*2;s_sub+sub*3];    
    %index for one sub-set.
end

%modified BSREM global parameters
%(1)Upper bound of reconstructed voxel intensity:BSREM II.
% load projMat_min
% for g=1:G
%     temp=sino(:,:,:,g)./pj;
%     int_U(g)=max(temp(:));
% end
% int_U=max(int_U);clear pj temp
int_U=5e5;
%(2)lower bound
epsilon=eps;
% nz_sino=find(sino>0);
% epsilon=min(sino(nz_sino));
% temp=ones(64,64,64);
% y=zeros(ns,nr,T,G);

%attenuation weights
% tic
if attn==1
    load weight128_120mn
    load weight128_attn1.mat
    wp_wgt=wp_attnwgt;clear wp_attnwgt;
else
    load weight128_120mn %no attenuation
end
% log_name=['wei_time_n' num2str(noisenum+1) '.mat'];
% wei_time=toc;save(log_name,'wei_time');
%load gbk64
% tic
load gbk128_beacon gb_temp

%(3)sub-seted denominator
pxn=128;tot_ang=120;
ang_span=42:109; %only 42~109 angular views nonzero
rsc=zeros(nx,ny,nz,'single');
for n_ind=1:T
    n=ang_span(n_ind);
    if n<tot_ang/2+1
        rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((n-1)*pxn+1:n*pxn),...
            wp_ipxl((n-1)*pxn+1:n*pxn),wp_wgt((n-1)*pxn+1:n*pxn),blur,gb_temp,tot_ang);
    else
        m=n-tot_ang/2;
        if attn==1%attenuation weights
            rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((m-1)*pxn+1:m*pxn),...
                wp_ipxl((m-1)*pxn+1:m*pxn),wp_wgt((n-1)*pxn+1:n*pxn),blur,gb_temp,tot_ang);
        else
            rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((m-1)*pxn+1:m*pxn),...
                wp_ipxl((m-1)*pxn+1:m*pxn),wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp,tot_ang);
        end
    end
end
zz_rsc=find(rsc==0);
rsc=rsc/sub;
%load motion compensation matrix: save time with the price of memory: 1GB
%extra memory needed. Save 16(g)*16(sub)*10(iter)*1.5(sec, load time)=3840 seconds

if gbeta~=0
    for g=1:G
        motionfile=['n4dMM' num2str(g) '_real.mat'];
        load(motionfile);
        M{g}=MM;
    end
    clear MM
end

load roi128
alpha0=1;
%objective function:time consuming. use for parameter determine.
%Afterwards not necessary! Necessary for adaptive alpha determine.
k=1;%Don't use it for real data OF; Revision not finished!!! May 02, 2007
if of_tag==1
    ps=zeros(G,1,'single');
    pg=zeros(G,1,'single');
    g_sub=zeros(ns,nr,T,G,'single');
    for g=1:G
        for m=1:T
            if m<tot_ang/2+1
                g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                    wp_vray((m-1)*pxn+1:m*pxn),...
                    wp_ipxl((m-1)*pxn+1:m*pxn),...
                    wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp,tot_ang);
            else
                comp_m=m-tot_ang/2;
                if attn==1%attenuation weights
                    g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                        wp_vray((comp_m-1)*pxn+1:comp_m*pxn),...
                        wp_ipxl((comp_m-1)*pxn+1:comp_m*pxn),...
                        wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp,tot_ang);
                else
                    g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                        wp_vray((comp_m-1)*pxn+1:comp_m*pxn),...
                        wp_ipxl((comp_m-1)*pxn+1:comp_m*pxn),...
                        wp_wgt((comp_m-1)*pxn+1:comp_m*pxn),blur,gb_temp,tot_ang);
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
st_ang=41;
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
                act_ang=st_ang+index(m,s);
                if act_ang<tot_ang/2+1
                    g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                        wp_vray((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_ipxl((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp,tot_ang);
                else
                    comp_ind=act_ang-tot_ang/2;
                    if attn==1%attenuation weights
                        g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp,tot_ang);
                    else
                        g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((comp_ind-1)*pxn+1:comp_ind*pxn),blur,gb_temp,tot_ang);
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
                if act_ang<tot_ang/2+1
                    temp=temp+back3d_saALL(g_sub,act_ang,...
                        wp_vray((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_ipxl((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp,tot_ang);
                else
                    comp_ind=act_ang-tot_ang/2;
                    if attn==1%attenuation weights
                        temp=temp+back3d_saALL(g_sub,act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp,tot_ang);
                    else
                        temp=temp+back3d_saALL(g_sub,act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((comp_ind-1)*pxn+1:comp_ind*pxn),blur,gb_temp,tot_ang);
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
                if m<tot_ang/2+1
                    g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                        wp_vray((m-1)*pxn+1:m*pxn),...
                        wp_ipxl((m-1)*pxn+1:m*pxn),...
                        wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp,tot_ang);
                else
                    comp_m=m-tot_ang/2;
                    if attn==1%attenuation weights
                        g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                            wp_vray((comp_m-1)*pxn+1:comp_m*pxn),...
                            wp_ipxl((comp_m-1)*pxn+1:comp_m*pxn),...
                            wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp,tot_ang);
                    else
                        g_sub(:,:,m,g)=proj3d_saALL(x(:,:,:,g),m,...
                            wp_vray((comp_m-1)*pxn+1:comp_m*pxn),...
                            wp_ipxl((comp_m-1)*pxn+1:comp_m*pxn),...
                            wp_wgt((comp_m-1)*pxn+1:comp_m*pxn),blur,gb_temp,tot_ang);
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