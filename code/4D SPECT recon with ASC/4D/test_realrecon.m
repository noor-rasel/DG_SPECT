%BSREM test
load real_gate1_sino temp
sino=temp(23:86,:,:);
x=ones(128,128,64);
iter=5;of_tag=0;sbeta=0;gbeta=0;blur=1;attn=0;scat=0;
tic,x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);toc
%11 min: 1 gate 5 iterations
save osem5 x

attn=1;

tic,x=mbsrem4dv3_real(sino,x,iter,of_tag,sbeta,gbeta,blur,attn,scat);toc
%12.5 min: 1 gate 5 iterations
save osem5_attn x
