function MM4dgen_noise30(n)
%from MM4dgen.m

%filename=['best spatial with AS correction noise #' num2str(n+1)];
%filename=['Im_ncatles_A_s5t0_n' num2str(n+1) '.mat']; for lesion
% filename=['./data_16g_5st0/Im_ncat16g_A_s4t0_n' num2str(n+1) '.mat']; %for normal
%filename=['./data_lesk_5st0/Im_ncatlesk_A_s4t0_n' num2str(n+1) '.mat']; %for King's LAD lesion

%filename=['./data_d64l_5st0/Im_d64l_A_s6t0_n' num2str(n+1) '.mat']; %256-->64; and Septal lesion
filename=['./data_d64_5st0/Im_d64_A_s4t0_n' num2str(n+1) '.mat']; %256-->64; and Septal lesion

load(filename)
%index in full  64^3
minx=24;maxx=48;
miny=16;maxy=39;
minz=29;maxz=48;
%index in sampled [23:52,16:43,29:48]
sminx=2;smaxx=26;
sminy=1;smaxy=24;
sminz=1;smaxz=20;
imx=zeros(64,64,64,16);
imx(minx:maxx,miny:maxy,minz:maxz,:)=Im_maps(sminx:smaxx,sminy:smaxy,sminz:smaxz,:);

G=16;
[x,y]=meshgrid(-3:3);FWHM=2;
gauss_kern=exp(-(x.^2+y.^2)*4*log(2)/FWHM^2);%if no "*4", it's HWHM.
gauss_kern=gauss_kern/sum(gauss_kern(:));%normalize
for g=1:G
    yin(:,:,:,g)=convn(imx(:,:,:,g),gauss_kern,'same');%(1.6207)%ones(3,3,3)/27
end

for i=1:G-1
    [vx,vy,vz]=motionele3d(yin(:,:,:,i:i+1),50,0.5);%alpha(n)
    vxx=zeros(size(vx));vyy=vxx;vzz=vxx;
    vxx(minx:maxx,miny:maxy,minz:maxz)=vx(minx:maxx,miny:maxy,minz:maxz);
    vyy(minx:maxx,miny:maxy,minz:maxz)=vy(minx:maxx,miny:maxy,minz:maxz);%vz-->vy Mar 14,2007
    vzz(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
    M1(:,:,:,:,i)=cat(4,vxx,vyy,vzz);
end
[vx,vy,vz]=motionele3d(yin(:,:,:,[G 1]),50,0.5);
vxx=zeros(size(vx));vyy=vxx;vzz=vxx;
vxx(minx:maxx,miny:maxy,minz:maxz)=vx(minx:maxx,miny:maxy,minz:maxz);
vyy(minx:maxx,miny:maxy,minz:maxz)=vy(minx:maxx,miny:maxy,minz:maxz);
vzz(minx:maxx,miny:maxy,minz:maxz)=vz(minx:maxx,miny:maxy,minz:maxz);
M1(:,:,:,:,G)=cat(4,vxx,vyy,vzz);
%save 4dNCAT_noiseM1 M1

%load 4dNCAT_noiseM1 M1
%motion for 5 gates
G=16;
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;%ind5(:,3)=[];
for gate=1:16
    MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
    filename=['n4dMM' num2str(gate) '_n' num2str(n)];
    save(filename,'MM');
end