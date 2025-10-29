function umdata_angle(s)

load real_gate1_sino temp
s_ang=1:53;
e_ang=68:120;
px_num=1;
rim=zeros(128,128,128);

rim=fbp_3dMALL(temp,px_num,s_ang(s+1),e_ang(s+1));

s_name=['fbp_realdata' num2str(s+1) '.mat'];
    
save(s_name,'rim')

%%%read MCAT
fid=fopen('mcat64_em_1','rb');
temp=fread(fid,'single');
fclose(fid);
temp=reshape(temp,[64 64 64]);
pha128=temp(:,:,35);
[xi,yi]=meshgrid(1:63/127:64,1:63/127:64);
pha128=interp2(pha128,xi,yi);

x=repmat(pha128,[1 1 128]);
y=zeros(128,128,60);
for n=1:60
    y(:,:,n)=proj3d_saALL(x,n,wp_vray((n-1)*128+1:n*128),wp_ipxl((n-1)*128+1:n*128),wp_wgt((n-1)*128+1:n*128),0,0,120);
end

rx=fbp_3dMALL(y,1,1,60);