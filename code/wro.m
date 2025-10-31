%extract weights from weigen3d.c output files
%reorgnized for proj3d_sa.m input
%circular orbit only
tic
fp=fopen('wgt_p64d64s64np64_3D','r');
ndet=64;nangle=64;
nrays=ndet*nangle;
%4096 rays; only 2048 needed because of symmetry
%projection on transaxial 2D slices! 
nrays_2=nrays/2;

%ray_header(1,nrays_2)=struct('pixels',[],'dcp',[],'hd_ptr',[]);%/* <--important structure */
ray_header_pixels=zeros(nrays_2,1);
ray_header_dcp=zeros(nrays_2,1);
wp_ipxl=cell(nrays_2,1);
wp_vray=cell(nrays_2,1);
wp_wgt=cell(nrays_2,1);k1=0;k2=0;k3=0;k4=0;k5=0;
for i=1:nrays_2
    npxls=fread(fp,1,'int');k1=k1+1;%4 bytes
    ray_header_pixels(i)=npxls;
    ray_header_dcp(i)=fread(fp,1,'single');k2=k2+1; %4 bytes/* distance from projection to first pixel */
    if (npxls ~= 0)
        wp=zeros(3,npxls);
        %wp(1,npxls)=struct('wgt',[],'ipxl',[],'vray',[]); %same structure as 'wgt_ptr'
        for j=1:npxls
            wp(1,j)=fread(fp,1,'uint16');k3=k3+1;%pixel index: 2 bytes
            wp(2,j)=fread(fp,1,'uint16');k4=k4+1;%ray index: 2 bytes
            wp(3,j)=fread(fp,1,'single');k5=k5+1;%weight: 4 bytes
        end
        %ray_header(i).hd_ptr=wp;
        wp_ipxl(i)={wp(1,:)};
        wp_vray(i)={wp(2,:)};
        wp_wgt(i)={wp(3,:)};
    end
end
fclose(fp);toc
save weight64 wp_vray wp_ipxl wp_wgt
%k1*4+k2*4+k3*2+k4*2+k5*4==1803584! Match size of 'wgt_p64d64s64np64_3D'.
%12sec, much smaller memory occupancy than 'struct'
%reorgnize for different angle weights in 2D transaxial slice!!!
% %pixel positions
% for n=1:64
%     wp=wp_ipxl{8*64+n};temp=zeros(64);temp(wp)=1;dsp(temp);pause;close;
% end
% %planes
% for n=1:64
%     wp_v=wp_vray{8*64+n};
%     if ~isempty(wp_v)
% %         wp_v_min(n)=wp_v(1);
% %         wp_v_max(n)=wp_v(end);
%     temp=zeros(64);
%     temp(wp_v(1)+1,:)=1;temp(wp_v(end)+1,:)=1;
%     dsp(temp);pause;close;
%     end
% end
% %weights
% wp=wp_ipxl{8*64+n};
% wp_w=wp_wgt{8*64+n};
% for n=1:152,temp(wp(n))=temp(wp(n))+wp_w(n);end
% wp_v(1:10)
% wp(1:10)
% wp_w(1:10)
% wp_v_max=zeros(2048,1);for n=1:2048,wp_v=wp_vray{n};if ~isempty(wp_v), wp_v_max(n)=max(wp_v);end;end
% %max plane # = 0-63(+1)
% %weight at the same (pixel) location is divided into different planes!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nplane=64;
% wp=zeros(64^2*nplane,nangle);
% for ang=1:nangle/2
%     for k=1:ndet
%         %to build first 32 planes (same distance to the detector) for each angle
%         wp_v=wp_vray{(ang-1)*64+k};%0-31
%         if ~isempty(wp_v)
%             %wp_v_max(ang,k)=max(wp_v);wp_v_min(ang,k)=min(wp_v);
%             wp_i=wp_ipxl{(ang-1)*64+k};dc=64*ones(size(wp_v));
%             iy=floor(wp_i./dc);%iy=iy+1;
%             ix=mod(wp_i,dc);ix=ix+1;
%             wp_i=iy.*dc+ix;
%             ind=wp_v*64^2+wp_i;
%             wp(ind,ang)=wp(ind,ang)+wp_wgt{(ang-1)*64+k}';
%         end
%     end
% end
% wp(wp<1e-6)=0;
% wp=reshape(wp,[64^2,64,64]);
% %(1) plane 1 have very small residual values (10^-7), set threshold 1e-6 to
% %get rid of them!
% %(2) some abnormal vaule at rightmost end, ignored!
% %(3) the next half projections are trivial since the circular orbit used (means symmetric!)
% %ray-detector information is still needed!!!
% for ang=nangle/2+1:nangle    
%     ang2=ang-nangle/2;
%     wp(:,:,ang)=wp(:,[64:-1:1],ang2);
% end
%Now we have all weights (wp 4096*64*64) for 64 angle positions 
%aligned on 64 planes at each angle. Applicable for all transaxial slices!

%cut down weights smaller than threshold
load weight64
threshold=0;%1e-4;
no=0;nm=0;
for j=1:2048
    wp_v=wp_vray{j};
    wp_i=wp_ipxl{j};
    wp_w=wp_wgt{j};
    if ~isempty(wp_w)
        L1=length(wp_v);no=no+L1;
        ind=find(wp_w>=threshold);
        wp_v=wp_v(ind);
        wp_i=wp_i(ind);
        wp_w=wp_w(ind);
        %modify plane:0-63 => 1-64
        wp_v=wp_v+1;
        %modify pixel locations: shift one row down!
        L=length(ind);nm=nm+L;
%         dc=64*ones(1,L);
%         iy=floor(wp_i/64);
%         ix=mod(wp_i,dc);ix=ix+1;
%         wp_i=iy*64+ix;
        wp_i=wp_i+1;
        wp_vray{j}=wp_v;
        wp_ipxl{j}=wp_i;
        wp_wgt{j}=wp_w;
    end
end% get rid of 2.78% weights!
%normalize!!! Mar. 23,2006
for j=1:32
    for i=1:64
        wp=wp_wgt{(j-1)*64+i};
        wp_S(i)=sum(wp);
    end
    wp_M(j)=max(wp_S);
end%find peak for each angle
for j=1:32
    for i=1:64
        wp_wgt{(j-1)*64+i}=wp_wgt{(j-1)*64+i}/max(wp_M);
    end
end%normalize
        
save weight64_mn wp_vray wp_ipxl wp_wgt
% for n=1:2048,wp=wp_wgt{n};gnn(n)=sum(wp);end
% k=0;
% for n=1:2048
%     wp=wp_ipxl{n};
%     if ~isempty(wp)
%         k=k+1;
%         min_ipxl(k)=min(wp);max_ipxl(k)=max(wp);
%     end
% end
% k=0;
% for n=1:2048
%     wp=wp_vray{n};
%     if ~isempty(wp)
%         k=k+1;
%         min_vray(k)=min(wp);max_vray(k)=max(wp);
%     end
% end
% %vray: 0-63
% fp=fopen('roi_p64d64s64np64_3D','r');
% roi=fread(fp,'uint16');
% fclose(fp);
% ind=find(roi>0);
% %min(ind)=90, min(min_ipxl)=89