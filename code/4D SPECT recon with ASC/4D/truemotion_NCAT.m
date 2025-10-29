%get true motion
filename='D:\imagereconstruction\4D\thomas_motion_ncat\ncat_mov16_vec_4_5.txt';
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
C = textscan(fid,'%*s %*s %f %f %f %*s %f %f %f %*s %f %f %f');
fclose(fid);
Xref = C{1};
Yref = C{2};
Zref = C{3};
Xtar = C{4};
Ytar = C{5};
Ztar = C{6};
Vtx = C{7};
Vty = C{8};
Vtz = C{9};
[lgr,ind]=size(Xref);
[x,y,z]=meshgrid(1:64);
load ncat4D16g_input ncat_phantom
temp=interp3(x,y,z,ncat_phantom(:,:,:,1),Xref,Yref,Zref);
Xmax=ceil(max(Xref));Xmin=floor(min(Xref));
Ymax=ceil(max(Yref));Ymin=floor(min(Yref));
Zmax=ceil(max(Zref));Zmin=floor(min(Zref));
dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1))
%how to use short motion vectors to compose sparse MC matrix?!
%ATTN: the reference frame also needs interpolation!!

%Easy solution to change reconstruction program as follows
%0.016*16*16*10*5/2=102.4 sec in cluster
%we only need to know each frames correspondece! very litter memory needed.
%typical 2079 components for each frame
betat=1;
for g=1:2,
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,[64 64 64]);
    ncat_phantom(:,:,:,g)=temp;
end
% Caution: Either transpose images, or switch Xref(row) and Yref(column) positions to
% get right interpolation results!
% for n=1:64
% 	ncat_phantom(:,:,n,1)=ncat_phantom(:,:,n,1)';
%     ncat_phantom(:,:,n,2)=ncat_phantom(:,:,n,2)';
% end
temp1=interp3(x,y,z,ncat_phantom(:,:,:,1),Yref,Xref,Zref);
temp2=interp3(x,y,z,ncat_phantom(:,:,:,2),Ytar,Xtar,Ztar);
temp3=interp3(x,y,z,ncat_phantom(:,:,:,2),Yref,Xref,Zref);
temp4=interp3(x,y,z,ncat_phantom(:,:,:,1),Yref,Xtar,Ztar);
prior=betat*(temp1-temp2)'*(temp1-temp2);%8.5823e+004:linear; 1.4421e+005: spline, very slow
prior=betat*(temp1-temp3)'*(temp1-temp3);%1.0244e+005
prior=betat*(temp3-temp4)'*(temp3-temp4);
for n=1:64
	ncat_phantom(:,:,n,1)=ncat_phantom(:,:,n,1)';
    ncat_phantom(:,:,n,2)=ncat_phantom(:,:,n,2)';
end
temp1=ncat_phantom(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1);dsp(temp1)
temp2=ncat_phantom(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,2);dsp(temp2)
prior=betat*(temp1(:)-temp2(:))'*(temp1(:)-temp2(:));%4.6098e+005
%2D display
% [slice,sZi]=sort(floor(Zref));
% ind=find(slice==35);
% figure,quiver(Xref(sZi(ind)),Yref(sZi(ind)),Vtx(sZi(ind)),Vty(sZi(ind)),0),axis ij
ind=find(Zref>=34.5&Zref<35.5);
figure,quiver(Xref(ind),Yref(ind),Vtx(ind),Vty(ind),0),axis ij
temp=ncat_phantom(:,:,35,1);
temp=mat2gray(temp);temp=1-temp/2;
figure,imagesc(temp',[0 1]),colormap(gray),axis equal,axis tight
hold on,quiver(Xref(ind)+1,Yref(ind),Vtx(ind),Vty(ind),0.6,'k'),axis([Xmin Xmax+1 Ymin Ymax]),axis off
hold on,plot(Xref(ind)+1,Yref(ind),'kx'),axis([Xmin Xmax+1 Ymin Ymax])

ind=find(Ztar>=35&Ztar<36);
temp=ncat_phantom(:,:,35,2);temp=mat2gray(temp);temp=1-temp/2;
figure,imagesc(temp',[0 1]),colormap(gray),axis equal,axis tight
hold on,plot(Xtar(ind)+1,Ytar(ind),'kx'),axis([Xmin Xmax+1 Ymin Ymax])

dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,1)')
dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,2)')
%3D
figure,quiver3(Xref,Yref,Zref,Vtx,Vty,Vtz,0)
%get each voxel
xr=round(Xref);
yr=round(Yref);
zr=round(Zref);figure,plot3(xr,yr,zr,'.')
I=sub2ind([64 64 64],yr,xr,zr);
[sI,ind]=sort(I);
k=0;xI_buff=0;nI_count=0;
for n=1:length(I)
    if sI(n)~=xI_buff
        k=k+1;
        xI_buff=sI(n);
        xI(k)=xI_buff;        
        nI_count(k)=1;                
    else
        nI_count(k)=nI_count(k)+1;        
    end
end
[X,Y,Z]=ind2sub([64 64 64],xI);
figure,plot3(X,Y,Z,'.')

%find vertices
m=64;n=64;k=64;
load roi
temp=ones(64,64,64);
for g=1:16
    if g~=16
        filename=['D:\imagereconstruction\4D\thomas_motion_ncat\ncat_mov16_vec_' num2str(g) '_' num2str(g+1) '.txt'];
    else
        filename=['D:\imagereconstruction\4D\thomas_motion_ncat\ncat_mov16_vec_' num2str(g) '_' num2str(1) '.txt'];
    end
    fid = fopen(filename);
    fgetl(fid);
    fgetl(fid);
    C = textscan(fid,'%*s %*s %f %f %f %*s %f %f %f %*s %f %f %f');
    fclose(fid);
    Xref{g} = C{1};
    Yref{g} = C{2};
    Zref{g} = C{3};
    xf=floor(Xref{g});yf=floor(Yref{g});zf=floor(Zref{g});
    xc=ceil(Xref{g});yc=ceil(Yref{g});zc=ceil(Zref{g});
    I=zeros(length(xf),8);
    I(:,1)=sub2ind([m n k],yf,xf,zf);
    I(:,2)=sub2ind([m n k],yf,xc,zf);
    I(:,3)=sub2ind([m n k],yc,xf,zf);
    I(:,4)=sub2ind([m n k],yc,xc,zf);
    I(:,5)=sub2ind([m n k],yf,xf,zc);
    I(:,6)=sub2ind([m n k],yf,xc,zc);
    I(:,7)=sub2ind([m n k],yc,xf,zc);
    I(:,8)=sub2ind([m n k],yc,xc,zc);
    for v=1:8
        temp(I(:,v))=0;
    end
end
temp=temp.*repmat(roi,[1 1 64]);
ind_out=find(temp>0);
for g=1:16
    len_heart(g)=length(Xref{g});
%     len_out(g)=length(ind{g});
end
save NCAT_motion_index Xref Yref Zref ind_out
%Thomas: not right?!
load ghost_ncat_mov16_vec_1_2.mat
for g=1:2,
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,[64 64 64]);
    ncat_phantom(:,:,:,g)=temp;
end
dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,1).*heart_ref(Xmin:Xmax,Ymin:Ymax,35));
dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,2).*heart_tar(Xmin:Xmax,Ymin:Ymax,35));
figure,quiver(Vx(Xmin:Xmax,Ymin:Ymax,35),Vy(Xmin:Xmax,Ymin:Ymax,35),0),axis ij
[x,y,z]=meshgrid(1:64);
xx=x+Vx;
yy=y+Vy;
zz=z+Vz;
temp1=ncat_phantom(:,:,:,2);
temp2=interp3(x,y,z,ncat_phantom(:,:,:,1),xx,yy,zz);
temp3=ncat_phantom(:,:,:,1);
betat=1;
prior=betat*(temp1(:)-temp2(:))'*(temp1(:)-temp2(:))
prior=betat*(temp1(:)-temp3(:))'*(temp1(:)-temp3(:))
%compare ideal OF and truth
load 4dNCAT_orgM1
M1=M1(:,:,:,:,1);
% [i,j]=meshgrid(Xmin:Xmax,Ymin:Ymax);
% figure,imagesc(ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,1)),colormap(gray),axis equal,axis tight
% hold on,quiver(M1(Xmin:Xmax,Ymin:Ymax,35,1),M1(Xmin:Xmax,Ymin:Ymax,35,2),0.6)%,axis ij
temp=ncat_phantom(Xmin:Xmax,Ymin:Ymax,35,1);
temp=mat2gray(temp);temp=1-temp/2;
figure,imagesc(temp',[0 1]),colormap(gray),axis equal,axis tight,axis off
hold on,quiver(M1(Xmin:Xmax,Ymin:Ymax,35,2)',M1(Xmin:Xmax,Ymin:Ymax,35,1)',0.6,'k')%,axis ij

temp1=ncat_phantom(:,:,:,1);
xx=x;yy=y;zz=z;
xx(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=x(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)+M1(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1);xx(xx>64)=64;xx(xx<1)=1;
yy(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=y(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)+M1(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1);yy(yy>64)=64;yy(yy<1)=1;
zz(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=z(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)+M1(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1);zz(zz>64)=64;zz(zz<1)=1;
temp2=interp3(x,y,z,ncat_phantom(:,:,:,2),xx,yy,zz);
temp3=ncat_phantom(:,:,:,2);
prior=betat*(temp1(:)-temp2(:))'*(temp1(:)-temp2(:))%4.0079e+005
prior=betat*(temp1(:)-temp3(:))'*(temp1(:)-temp3(:))
%see true_MM4dgen.m and true_ncatmotion_batch.m 