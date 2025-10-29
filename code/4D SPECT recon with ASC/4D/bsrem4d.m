function [x,obj_func]=bsrem4d(sino,x,sub,iter,of_tag,sbeta,gbeta,blur)
%Original block sequential regularized EM.
%Slower, but better convergence than modified version.
%3D and 4D with motion compensation.
%Mingwu Jin, July 2006

if nargin<7, blur=1;end%default SPECT
[nx,ny,nz,G]=size(x);%G gates.
[ns,nr,T,G]=size(sino);%nr=nx=ny(bin 64); ns=nz(slice 64); T (angle: 64).

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

%BSREM global parameters

load weight64_mn
load gbk64


load roi
%main loop
for n=1:iter
    alpha=1/n;%0.015/n^.33*2;
    %most important for convergence! Need to adjust for different case.
    %1/(n/15+1) for 8 subsets in Ahn paper, good for EM 16 subsets 10 iter. 
    %1/(n*4+1) for sbeta=0.01. 5 iterations
    
    %objective function:time consuming. use for parameter determine.
    %Afterwards not necessary!
    if of_tag==1
        ps=zeros(G,1);
        pg=zeros(G,1);
        for g=1:G
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
                %prior 1: (Xi-1/Ni*sum(Xj))^2, Xj is the neighbor f Xi
                temp=imfilter(x(:,:,:,g),sfilt,'replicate','same');
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
        obj_func(n)=sum(sino(nz_g).*log(g_sub(nz_g))-g_sub(nz_g))-sum(ps)-sum(pg);
    else
        obj_func=[];
    end
    
    for s=1:sub
        y=zeros(nx,ny,nz,G);
        for g=1:G
            temp=zeros(nx,ny,nz);
            for m=1:sub_com
                sino_sub=zeros(nz,nx);
                g_sub=zeros(nz,nx);
                sino_sub=sino(:,:,index(m,s),g);
                if index(m,s)<33
                    g_sub=proj3d_sa(x(:,:,:,g),index(m,s),...
                        wp_vray((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_ipxl((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                else
                    comp_ind=index(m,s)-32;
                    g_sub=proj3d_sa(x(:,:,:,g),index(m,s),...
                        wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                        wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                        wp_wgt((comp_ind-1)*64+1:comp_ind*64),blur,gb_temp);
                end
                nnz=find(g_sub>eps);
                g_sub(nnz)=sino_sub(nnz)./g_sub(nnz)-1;                
                %derivative of objective function:backprojection of ratio
                if index(m,s)<33
                    temp=temp+back3d_sa(g_sub,index(m,s),...
                        wp_vray((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_ipxl((index(m,s)-1)*64+1:index(m,s)*64),...
                        wp_wgt((index(m,s)-1)*64+1:index(m,s)*64),blur,gb_temp);
                else
                    comp_ind=index(m,s)-32;
                    temp=temp+back3d_sa(g_sub,index(m,s),...
                        wp_vray((comp_ind-1)*64+1:comp_ind*64),...
                        wp_ipxl((comp_ind-1)*64+1:comp_ind*64),...
                        wp_wgt((comp_ind-1)*64+1:comp_ind*64),blur,gb_temp);
                end
            end
            %derivative of objective function:backprojection of ratio
            y(:,:,:,g)=temp;
        end
        %update
        x=x+alpha*x.*y;        
        if s==sub
            y=zeros(nx,ny,nz,G);
            for g=1:G
                if sbeta~=0
                    %derivative of objective function:spatial prior
                    %prior 1:
                    temp=imfilter(x(:,:,:,g),sfilt,'replicate','same');
                    temp=2*sbeta*imfilter(temp,sfilt,'replicate','same');
                    %prior 2:
                    %temp=2*sbeta*imfilter(x(:,:,:,g),sfilt,'replicate','same')/6;
                    y(:,:,:,g)=y(:,:,:,g)-temp;
                end
                if gbeta~=0
                    %derivative of objective function:temporal prior (5 gates)
                    temp=x(:,:,:,ind5(g,:));
                    motionfile=['MM4dgate' num2str(g)];
                    load(motionfile);
                    temp=gbeta*MM*temp(:);%<=> x-MM*x
                    temp=reshape(temp,[nx ny nz]);
                    y(:,:,:,g)=y(:,:,:,g)-temp;
                end
            end
            x=x+alpha*x.*y;
            x(x<=0)=eps;
        end
    end
end