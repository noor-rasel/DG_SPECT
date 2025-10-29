function test_realASC_optbetas(s)

betas=[0 1e-5 5e-5 1e-4 3e-4];
% betas=[1e-4 3e-4];
%1e-3 too big; 1e-6 too small
[k,g]=ind2sub([5 8],s+1);

load real_16g_sino
sino=sino(:,:,:,g);tc=sum(sino(:));
x=ones(128,128,64);
iter=10;of_tag=0;sbeta=betas(k);gbeta=0;blur=1;attn=1;scat=0;
load real_16g_scat scat
scat=scat(:,:,:,g);

x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);
x=x*tc/sum(x(:));

filename=['realg' num2str(g) '_mapASC_s' num2str(k) '.mat'];
save(filename,'x');