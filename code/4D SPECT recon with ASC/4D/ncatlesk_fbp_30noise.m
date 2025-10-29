function ncatlesk_fbp_30noise(k)
%16 gates NCAT phantom FBP

filename=['ncatlesk_n' num2str(k+1) '.mat'];
load(filename)

Im_rawfbp=fbp_3dM(sinoT);

Im_fbp=clinicfilt3d_yyjin(Im_rawfbp,2.4,.4);
Im_fbp=Im_fbp*8e6/sum(Im_fbp(:));
Im_fbp=Im_fbp(23:52,16:43,29:48,:);

filename=['./data_lesk_fbp/Im_fbp_n' num2str(k+1) '.mat'];
save(filename,'Im_fbp');
