function test_real_fbp(s)

g=s+1;

load real_16g_sino
sino=sino(:,:,:,g);tc=sum(sino(:));
x=ones(128,128,64);px_num=1;

x=fbp_3dMALL(sino,px_num,42,42+67);
x=clinicfilt3d_yyjin(x,2.4,.4);
x=x*tc/sum(x(:));

s_name=['realg' num2str(g) '_fbp.mat'];
    
save(s_name,'x')
