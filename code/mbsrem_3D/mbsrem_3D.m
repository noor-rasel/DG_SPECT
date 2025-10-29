function [x]=mbsrem_3D(sino,x,sub,iter,sbeta,blur,attn,scat,acq_time,pid)
%Modified block sequential regularized EM.
%3D and 4D with motion compensation.
%Mingwu Jin, July 2006
%add attenuation correction (AC) and scatter correction (SC) Oct 07,2006.
%of_tag=1: used to calcualte MAP objective function and guarantee its
%increasing.
%single precision: Nov. 09, 2006
%add true motion: Nov. 13, 2006
%dedicated for IRIX clincal data
%using different ror for each projection: Aug, 2013
%imcorporate the acquisition time variance for different angle.


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
 
% for a spatial blur.
 mask = ones(3,3,3);
 mask(2,2,2) = 3;
 mask = mask./sum(mask(:));

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


%generate subset index:"index";
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
if attn==1
    load weight128_mn
    load (['/Users/wenyuanqi/Desktop/clinical/shared_data/reference/' pid '/' pid '_weight128_attn.mat']);
    wp_wgt=wp_attnwgt;clear wp_attnwgt;
else
    load weight128_mn %no attenuation
end
% log_name=['wei_time_n' num2str(noisenum+1) '.mat'];
% wei_time=toc;save(log_name,'wei_time');
%load gbk64
load (['/Users/wenyuanqi/Desktop/clinical/shared_data/reference/' pid '/' pid '_gbk128.mat']);

%(3)sub-seted denominator
pxn=128;tot_ang=128;
ang_span=49:112; 
rsc=zeros(nx,ny,nz,'single');
for n_ind=1:T
    n=ang_span(n_ind);
    
    clear gb_temp0;
    for ppi = 1:128
       gb_temp0{ppi} = gb_temp{ppi,n_ind}; 
    end
    
    if n<tot_ang/2+1
        rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((n-1)*pxn+1:n*pxn),...
            wp_ipxl((n-1)*pxn+1:n*pxn),wp_wgt((n-1)*pxn+1:n*pxn),blur,gb_temp0,tot_ang);
    else
        m=n-tot_ang/2;
        if attn==1%attenuation weights
            rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((m-1)*pxn+1:m*pxn),...
                wp_ipxl((m-1)*pxn+1:m*pxn),wp_wgt((n-1)*pxn+1:n*pxn),blur,gb_temp0,tot_ang);
        else
            rsc=rsc+back3d_saALL(ones(nz,nx,'single'),n,wp_vray((m-1)*pxn+1:m*pxn),...
                wp_ipxl((m-1)*pxn+1:m*pxn),wp_wgt((m-1)*pxn+1:m*pxn),blur,gb_temp0,tot_ang);
        end
    end
end
rsc=rsc/sub;
%load motion compensation matrix: save time with the price of memory: 1GB
%extra memory needed. Save 16(g)*16(sub)*10(iter)*1.5(sec, load time)=3840 seconds

load roi128
alpha0=1;
%objective function:time consuming. use for parameter determine.
%Afterwards not necessary! Necessary for adaptive alpha determine.
k=1;%Don't use it for real data OF; Revision not finished!!! May 02, 2007

%main loop
st_ang=48;

for n=1:iter
    alpha=alpha0/(k+1);%0.015/n^.33*2;
    %most important for convergence! Need to adjust for different case.
    %1/(n/15+1) for 8 subsets in Ahn paper, good for EM 16 subsets 10 iter.
    %1/(n*4+1) for sbeta=0.01. 5 iterations
    for s=1:sub
        y=zeros(nx,ny,nz,G);
        for g=1:G
            temp=zeros(nx,ny,nz);
            for m=1:sub_com
            if acq_time(index(m,s),g)>0.0001
                                  
                g_sub=zeros(nz,nx);
                sino_sub=sino(:,:,index(m,s),g);
                act_ang=st_ang+index(m,s);
                
                  clear gb_temp0;
                    for ppi = 1:128
                       gb_temp0{ppi} = gb_temp{ppi,index(m,s)}; 
                    end
     
                
                if act_ang<tot_ang/2+1
                    g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                        wp_vray((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_ipxl((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp0,tot_ang);
                else
                    comp_ind=act_ang-tot_ang/2;
                    if attn==1%attenuation weights
                        g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp0,tot_ang);
                    else
                        g_sub=proj3d_saALL(x(:,:,:,g),act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((comp_ind-1)*pxn+1:comp_ind*pxn),blur,gb_temp0,tot_ang);
                    end
                end
                %scatter compensation
                if length(scat)>1
                    g_sub=g_sub+scat(:,:,index(m,s),g)/acq_time(index(m,s),g);
                end
                nnz=find(g_sub>epsilon);
          
                g_sub(nnz)=sino_sub(nnz)./g_sub(nnz)-acq_time(index(m,s),g);
                %derivative of objective function:backprojection of ratio
                if act_ang<tot_ang/2+1
                    temp=temp+back3d_saALL(g_sub,act_ang,...
                        wp_vray((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_ipxl((act_ang-1)*pxn+1:act_ang*pxn),...
                        wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp0,tot_ang);
                else
                    comp_ind=act_ang-tot_ang/2;
                    if attn==1%attenuation weights
                        temp=temp+back3d_saALL(g_sub,act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((act_ang-1)*pxn+1:act_ang*pxn),blur,gb_temp0,tot_ang);
                    else
                        temp=temp+back3d_saALL(g_sub,act_ang,...
                            wp_vray((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_ipxl((comp_ind-1)*pxn+1:comp_ind*pxn),...
                            wp_wgt((comp_ind-1)*pxn+1:comp_ind*pxn),blur,gb_temp0,tot_ang);
                    end
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
          
            %step size of derivative
            temp=x(:,:,:,g);
            temp1=zeros(nx,ny,nz);
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
    k=k+1;
end