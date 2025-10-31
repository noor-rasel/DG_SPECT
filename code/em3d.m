function x=em3d(d,x,iter,blur)%,li,sbeta,gbeta,motion)

if nargin<4, blur=1;end%default PET
[nx,ny,nz,T,G]=size(x);
fp=fopen('roi_p64d64s64np64_3D','r');
roi=fread(fp,'uint16');
fclose(fp);
roi=repmat(roi,1,64);
ind=find(roi>0);clear roi%4096->3072
NC=0.01575;%? normalization need for error control!!!

rsc=ones(64,64,64);
rsc=back3d(rsc,blur,1,'-p',64,'-d',64,'-s',64,'-ror',28.5,'-np',64);
rsc=rsc*NC;

for k=1:iter
    %projection
    y=proj3d(x,blur,1,'-p',64,'-d',64,'-s',64,'-ror',28.5,'-np',64);y=y*NC;
    pro_nz=find(y>0);
    y(pro_nz)=d(pro_nz)./y(pro_nz);
    %backprojection
    t=back3d(y,blur,1,'-p',64,'-d',64,'-s',64,'-ror',28.5,'-np',64);t=t*NC;
    t(ind)=t(ind)./rsc(ind);
    x=x.*t;
end