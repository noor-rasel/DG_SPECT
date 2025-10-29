% Erwan Gravier
% Research Project
% 2 Oct. 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help file for weigen3d (C program)
% This program creates the necessary files for the forward (proj3d)
% and backward projection (back3d) to work.
% The original C program by R. L. Siddon has been modified to adapt it to the matlab
% environment.
%
% weigen3d('-p',pxls[64],'-d',dets[64],'-s',steps[64],'-np',planes[64],
%      '-rpx',[reconstruction pixel size, 0.712cm for 64*64 objects],
%      '-apx',[=rpx, acquisition pixel size, 0.712cm for 64*64 objects])
%
%           pxls: number of pixels in each direction.
%           dets: number of detectors in lateral and axial directions.
%           steps: number of projection angles (even number).
%           planes: number of planes in projections.
%           rpx: reconstruction pixel size.
%           apx: acquisition pixel size.
%
% No matlab output generated.
% Several files are generated:
%           - .rec_p64d64s64np64_3D (1KB)
%           - rec_p64d64s64np64_3D  (1KB)
%           - roi_p64d64s64np64_3D  (8KB)
%           - wgt_p64d64s64np64_3D  (1,762KB) weight file
%
% Attention!!! To get more information, type: weigen3d('help')
% Ex. weigen3d('-p',64,'-d',64,'-s',64,'-np',64,'-rpx',0.712,'-apx',0.712)
% Time required: 0.7660s on P4 2GHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%