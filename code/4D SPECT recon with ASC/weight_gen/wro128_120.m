tic
fp=fopen('wgt_p128d128s120np128_3D','r');
ndet=128;nangle=120;
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
save weight128_120 wp_vray wp_ipxl wp_wgt

load weight128_120
threshold=0;%1e-4;
no=0;nm=0;
for j=1:nrays_2
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
for j=1:60
    for i=1:128
        wp=wp_wgt{(j-1)*128+i};
        wp_S(i)=sum(wp);
    end
    wp_M(j)=max(wp_S);
end%find peak for each angle
for j=1:60
    for i=1:128
        wp_wgt{(j-1)*128+i}=wp_wgt{(j-1)*128+i}/max(wp_M);
    end
end%normalize
        
save weight128_120mn wp_vray wp_ipxl wp_wgt