function MM4dgen_real(n)
%from MM4dgen.m

filename=['real_data_mapsAC_oim.mat']; %for King's LAD lesion
load(filename)
%index in full  64^3
G=8;
miny=48;maxy=73;
minx=44;maxx=69;
minz=17;maxz=31;
%index in sampled [23:52,16:43,29:48]
yin=zeros(128,128,64,G);
yin(minx:maxx,miny:maxy,minz:maxz,:)=y(minx:maxx,miny:maxy,minz:maxz,:);


% for g=1:G
%     yin(:,:,:,g)=convn(imx(:,:,:,g),ones(3,3,3)/27,'same');%(1.6207)
% end

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
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;%ind5(:,3)=[];
for gate=1:G
    MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
    filename=['n4dMM' num2str(gate) '_real'];
    save(filename,'MM');
end