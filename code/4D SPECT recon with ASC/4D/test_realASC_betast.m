function test_realASC_betast(s)

betas=[0 1e-5 3e-5 5e-5 1e-4 3e-4];
betat=[1e-5 2e-5 3e-5 4e-5 5e-5 1e-4];
% betas=[1e-4 3e-4];
%1e-3 too big; 1e-6 too small
[k,t]=ind2sub([6 6],s+1);

load real_16g_sino
tc=sum(sino(:));
x=ones(128,128,64,8);
iter=10;of_tag=0;sbeta=betas(k);gbeta=betat(t);blur=1;attn=1;scat=0;
load real_16g_scat scat

x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);
x=x*tc/sum(x(:));

filename=['real8g_mapASC_s' num2str(k) 't' num2str(t) '.mat'];
save(filename,'x');