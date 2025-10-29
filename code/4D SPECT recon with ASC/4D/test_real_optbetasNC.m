function test_real_optbetasNC(s)

%betas=[0 1e-4 5e-4 1e-3 3e-3];
betas=3e-3;
g=s+1;

load real_16g_sino
sino=sino(:,:,:,g);tc=sum(sino(:));
x=ones(128,128,64);
iter=10;of_tag=0;sbeta=betas;gbeta=0;blur=1;attn=0;scat=0;

x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);
x=x*tc/sum(x(:));

filename=['realg' num2str(g) '_mapNC_s' num2str(5) '.mat'];
save(filename,'x');