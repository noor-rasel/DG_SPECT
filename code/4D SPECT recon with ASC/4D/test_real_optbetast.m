function test_real_optbetast(s)


load real_16g_sino
% sino=sino(:,:,:,g);
betas=[1e-5 3e-5 5e-5];
betat=[1e-5 2e-5 3e-5 4e-5 5e-5];
[k,t]=ind2sub([3 5],s+1);%14
tc=sum(sino(:));
x=ones(128,128,64,8);
iter=10;of_tag=0;sbeta=betas(k);gbeta=betat(t);blur=1;attn=1;scat=0;

x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);
x=x*tc/sum(x(:));

filename=['real8gate_map_s' num2str(k) 't' num2str(t) '.mat'];
save(filename,'x');
%real8gate_map_ST.mat%s: 5e-5; t: 5e-5