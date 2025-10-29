function MM4dgen_noise(gate)

load 4dNCAT_noiseM1 M1
%motion for 5 gates
G=16;
for g=1:G
    ind5(g,:)=mod((g-2:g+2),G);
end
ind5(ind5==0)=G;%ind5(:,3)=[];
gate=gate+1;
    MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
    filename=['n4dMM' num2str(gate)];
    save(filename,'MM');