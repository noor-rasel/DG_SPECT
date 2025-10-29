function MM4dgen(Im_maps)
%generate motion prediction matrix for each gate from adjacent 5 gates

G=8;
minx=24;maxx=48;
miny=16;maxy=39;
minz=29;maxz=48;
%index in sampled [23:52,16:43,29:48]
sminx=2;smaxx=26;
sminy=1;smaxy=24;
sminz=1;smaxz=20;
imx=zeros(64,64,64,G);
imx(minx:maxx,miny:maxy,minz:maxz,:)=Im_maps(sminx:smaxx,sminy:smaxy,sminz:smaxz,:);


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

%motion for 5 gates
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;
for gate=1:G
    MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
    filename=['n4dMM' num2str(gate)];
    save(filename,'MM');
end