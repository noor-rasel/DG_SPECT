%true NCAT motion generation and verification

minx=255;
miny=255;
minz=255;
maxx=0;
maxy=0;
maxz=0;
for fr=1:16
    [Vtx,Vty,Vtz,Xref,Yref,Zref]=load_motion_frame(fr);
    minx=floor(min(min(Xref)-1,minx));
    miny=floor(min(min(Yref)-1,miny));
    minz=floor(min(min(Zref)-1,minz));

    maxx=ceil(max(max(Xref)+1,maxx));
    maxy=ceil(max(max(Yref)+1,maxy));
    maxz=ceil(max(max(Zref)+1,maxz));

end

for fr=1:16;
    [Vtx,Vty,Vtz,Xref,Yref,Zref,Loc]=load_motion_frame(fr);
    fname=sprintf('myoncat16g_act_%i.bin',fr);
    fid=fopen(fname,'rb');
    F=fread(fid,'float');
    mask=reshape(F,64,64,64);
    mask=permute(mask,[ 2 1 3]);

    mask=convn(mask,ones(2,2,2),'same');

    fclose(fid);

    [Reg_x,Reg_y,Reg_z]=meshgrid(minx:maxx,miny:maxy,minz:maxz);

    Reg_x=Reg_x.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
    Reg_y=Reg_y.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
    Reg_z=Reg_z.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
    %
    Reg_x=Reg_x(Reg_x>0);
    Reg_y=Reg_y(Reg_y>0);
    Reg_z=Reg_z(Reg_z>0);

    [Reg_Vtx]=griddata3(Xref,Yref,Zref,Vtx,Reg_x,Reg_y,Reg_z,'nearest',[]);
    [Reg_Vty]=griddata3(Xref,Yref,Zref,Vty,Reg_x,Reg_y,Reg_z,'nearest',[]);
    [Reg_Vtz]=griddata3(Xref,Yref,Zref,Vtz,Reg_x,Reg_y,Reg_z,'nearest',[]);

    Reg_Vtx(isnan(Reg_Vtx))=0;
    Reg_Vty(isnan(Reg_Vty))=0;
    Reg_Vtz(isnan(Reg_Vtz))=0;
    %show motion
    Slice=36
    dspm(mask(:,:,Slice)')
    hold on
    IndPlot=((Zref>Slice-0.5)&(Zref<Slice+0.5));
    quiver(Xref(IndPlot),Yref(IndPlot),Vtx(IndPlot),Vty(IndPlot),0,'r')
    hold on
    Slice=Slice;
    IndPlot2=((Reg_z>=Slice)&(Reg_z<=Slice));
    quiver(Reg_x(IndPlot2),Reg_y(IndPlot2),Reg_Vtx(IndPlot2),Reg_Vty(IndPlot2),0,'y')
    hold off
    axis([minx,maxx,miny,maxy]);
eval(sprintf('save motion_frame_%i Reg_x Reg_y Reg_z Reg_Vtx Reg_Vty Reg_Vtz',fr))
end

%direct generate MC matrix
G=16;
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;ind5(:,3)=[];
for fr=1:G
    %current frame coordinates
    [Vtx,Vty,Vtz,Xref,Yref,Zref,Loc]=load_motion_frame(fr);
    %heart mask
    fname=sprintf('myoncat16g_act_%i.bin',fr);
    fid=fopen(fname,'rb');
    F=fread(fid,'float');
    mask=reshape(F,64,64,64);
    mask=permute(mask,[ 2 1 3]);
    mask=convn(mask,ones(2,2,2),'same');
    fclose(fid);

    [Reg_x,Reg_y,Reg_z]=meshgrid(minx:maxx,miny:maxy,minz:maxz);
%     Reg_ind=sub2ind([64 64 64],Reg_y,Reg_x,Reg_z);
    Reg_x=Reg_x.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
    Reg_y=Reg_y.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
    Reg_z=Reg_z.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);        
    MM=cell(4);
    for n=1:4        
        %target frame coordinates
        [Vtx,Vty,Vtz,Xtar,Ytar,Ztar,Loc]=load_motion_frame(ind5(fr,n));
        Vtx=Xtar-Xref;
        Vty=Ytar-Yref;
        Vtz=Ztar-Zref;
        [Reg_Vtx]=griddata3(Xref,Yref,Zref,Vtx,Reg_x,Reg_y,Reg_z,'nearest',{'Qt','Qbb','Qc'});
        [Reg_Vty]=griddata3(Xref,Yref,Zref,Vty,Reg_x,Reg_y,Reg_z,'nearest',{'Qt','Qbb','Qc'});
        [Reg_Vtz]=griddata3(Xref,Yref,Zref,Vtz,Reg_x,Reg_y,Reg_z,'nearest',{'Qt','Qbb','Qc'});
        Reg_Vtx(isnan(Reg_Vtx))=0;
        Reg_Vty(isnan(Reg_Vty))=0;
        Reg_Vtz(isnan(Reg_Vtz))=0;
        Reg_Vtx=Reg_Vtx.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
        Reg_Vty=Reg_Vty.*(mask(miny:maxy,minx:maxx,minz:maxz)>0);
        Reg_Vtz=Reg_Vtz.*(mask(miny:maxy,minx:maxx,minz:maxz)>0); 
        vx=zeros(64,64,64);vy=vx;vz=vx;
        vx(miny:maxy,minx:maxx,minz:maxz)=-Reg_Vtx;
        vy(miny:maxy,minx:maxx,minz:maxz)=-Reg_Vty;
        vz(miny:maxy,minx:maxx,minz:maxz)=-Reg_Vtz;
        tic,M1=motionmatrix3d(permute(vy,[2 1 3]),permute(vx,[2 1 3]),permute(vz,[2 1 3]));toc%161sec in P4
    end
    MM{n}=M1;
    MM=cat(2,MM{:});
    filename=['true_4dMM' num2str(g)];
    save(filename,'MM');
end
    
    
%verification
%single
fname=sprintf('myoncat16g_act_%i.bin',fr);
fid=fopen(fname,'rb');
cf=fread(fid,'float');
fclose(fid);
fname=sprintf('myoncat16g_act_%i.bin',ind5(fr,n));
fid=fopen(fname,'rb');
tf=fread(fid,'float');
fclose(fid);
mean((cf(:)-tf(:)).^2)
mean((cf(:)-(tf(:)'*M1)').^2)
cf=reshape(cf,64,64,64);
tf=reshape(tf,64,64,64);
cf=permute(cf,[2 1 3]);
tf=permute(tf,[2 1 3]);
% for sli=1:64,cf(:,:,sli)=cf(:,:,sli)';tf(:,:,sli)=tf(:,:,sli)';end
mcf=zeros(size(cf));
mcf(:)=(tf(:)'*M1)';
%multiple: 4
load true_4dMM1
for fr=1:16
    fname=sprintf('lesion16g_act_%i.bin',fr);
    fid=fopen(fname,'rb');
    cf=fread(fid,'float');
    fclose(fid);
    temp(:,fr)=cf;
end
for fr=1:16
    load(['true_4dMM' num2str(fr) '.mat']);
    temp1=temp(:,ind5(fr,:));
    mse_m(fr)=mean((temp(:,fr)-mean(temp1,2)).^2);
    mse_mc(fr)=mean((temp(:,fr)-MM*temp1(:)).^2);
end
%not very good, but still better than average!
%Jovan
% [Int_x,Int_y,Int_z]=meshgrid(1:64,1:64,1:64);
% int_Vtx=zeros(64,64,64);
% int_Vty=zeros(64,64,64);
% int_Vtz=zeros(64,64,64);
% int_Vtx(motion_mask(:)~=0)=Reg_Vtx;
% int_Vty(motion_mask(:)~=0)=Reg_Vty;
% int_Vtz(motion_mask(:)~=0)=Reg_Vtz;
% 
% 
% fim=interp3(Int_x,Int_y,Int_z,im,Int_y+int_Vty,Int_x+int_Vtx,Int_z+int_Vtz);
