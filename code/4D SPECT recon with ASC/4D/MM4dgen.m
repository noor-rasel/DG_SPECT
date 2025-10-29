function MM4dgen(gate)
% load ncat4D16g_input ncat_phantom
% x=ncat_phantom;clear ncat_phantom
% yin=zeros(size(x));
G=16;
% for g=1:G
%     yin(:,:,:,g)=convn(x(:,:,:,g),ones(3,3,3)/27,'same');
% end
% M1=zeros(64,64,64,3,G);
% for i=1:G-1
%     [vx,vy,vz]=motionele3d(yin(:,:,:,i:i+1),50,.1);
%     M1(:,:,:,:,i)=cat(4,vx,vy,vz);
% end
% [vx,vy,vz]=motionele3d(yin(:,:,:,[G 1]),50,.1);
% M1(:,:,:,:,G)=cat(4,vx,vy,vz);
% save 4dNCAT_orgM1 M1
load 4dNCAT_orgM1 M1
%motion for 5 gates
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;
gate=gate+8;
%for g=8:G
    MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
    filename=['4dMM' num2str(gate)];
    save(filename,'MM');
%end
%validate
% load ncat4D16g_input ncat_phantom
% x=ncat_phantom;clear ncat_phantom
% yin=zeros(size(x));
% G=16;
% for g=1:G
%     yin(:,:,:,g)=convn(x(:,:,:,g),ones(3,3,3)/27,'same');
% end
% clear x
% for g=1:G
%     ind5(g,:)=mod((g-2:g+2),G);
% end
% ind5(ind5==0)=G;
% load 4dMM1
%%(1) gamma=1, linear slop
% temp=yin(:,:,:,ind5(1,:));
% temp1=temp(:,:,:,3);
% ssd_M=sum((temp1(:)-MM*temp(:)).^2)%486.6536
%%(2) simple averge
% % MM0=cat(2,speye(64^3,64^3)/4,speye(64^3,64^3)/4);
% % MM0=cat(2,MM0,spalloc(64^3,64^3,0));
% % MM0=cat(2,MM0,speye(64^3,64^3)/4);
% % MM0=cat(2,MM0,speye(64^3,64^3)/4);
% % ssd_A=sum((temp1(:)-MM0*temp(:)).^2)%2.4300e+03
% tempA=mean(yin(:,:,:,[15 16 2 3]),4);
% ssd_A=sum((temp1(:)-tempA(:)).^2)
%%(3) gamma=0, same weights
% load 4dMM1
% temp=yin(:,:,:,ind5(1,:));
% G=5;g=3;gamma=1;
% k=1:G;%% conlusion: motion are useful!!!
% wei=abs(1-2*abs(g-k)/G).^gamma;
% wei(g)=0;wei=wei/sum(wei);
% for n=1:G
%     if n~=3
%         temp(:,:,:,n)=temp(:,:,:,n)/wei(n);
%     end
% end
% temp1=temp(:,:,:,3);
% ssd_M=sum((temp1(:)-MM*temp(:)/4).^2)%511.9315
%% conlusion: motion are useful!!! gamma=1 is best

%%%%%%%%%%%%%%%
%noisy motion
load ncatlestew_maps_ACopt.mat
minx=24;maxx=48;
miny=16;maxy=39;
minz=26;maxz=47;
imx=zeros(size(Im_maps));
imx(minx:maxx,miny:maxy,minz:maxz,:)=Im_maps(minx:maxx,miny:maxy,minz:maxz,:);
% yin=clinicfilt3d_yyjin(imx,3,0.2);% 0.2(1.6302) 0.3(1.6268) 0.4(1.6244)
% 0.6(1.6233); no filter(1.6218)
G=16;
for g=1:G
    yin(:,:,:,g)=convn(imx(:,:,:,g),ones(3,3,3)/27,'same');%(1.6207)
end
% for fr=1:16
%     fname=sprintf('./motion_frame/ncat16gles_act_%i.bin',fr);
%     fid=fopen(fname,'rb');
%     cf=fread(fid,'float');
%     fclose(fid);
%     temp(:,fr)=cf;
% end
% temp=reshape(temp,64,64,64,16);
% M1=zeros(64,64,64,3,G);
% [x,y,z]=meshgrid(1:64,1:64,1:64);
% alpha=0.1:0.1:1;
for i=1:G-1
%     cf=temp(:,:,:,i);tf=temp(:,:,:,i+1);
%     mse_a(i)=mean((cf(:)-tf(:)).^2);
% for n=1:10
    [vx,vy,vz]=motionele3d(yin(:,:,:,i:i+1),50,0.5);%alpha(n)
    vxx=zeros(size(vx));vyy=vxx;vzz=vxx;
    vxx(minx:maxx,miny:maxy,minz:maxz)=vx(minx:maxx,miny:maxy,minz:maxz);
    vyy(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
    vzz(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
   M1(:,:,:,:,i)=cat(4,vxx,vyy,vzz); 
%     tf=temp(:,:,:,i+1);
%     tf=interp3(x,y,z,tf,x+vxx,y+vyy,z+vzz);
%     mse_mcn(n)=mean((cf(:)-tf(:)).^2);
% end
end
[vx,vy,vz]=motionele3d(yin(:,:,:,[G 1]),50,0.5);
    vxx=zeros(size(vx));vyy=vxx;vzz=vxx;
vxx(minx:maxx,miny:maxy,minz:maxz)=vx(minx:maxx,miny:maxy,minz:maxz);
    vyy(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
    vzz(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
M1(:,:,:,:,G)=cat(4,vxx,vyy,vzz);
save 4dNCAT_noiseM1 M1