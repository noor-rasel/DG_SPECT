function [Vtx,Vty,Vtz,Xref,Yref,Zref,Loc]=load_motion_frame(fr)
filename=['D:\imagereconstruction\4D\thomas_motion_ncat\ncat_mov16_vec_'...
    num2str(fr) '_' num2str(mod(fr,16)+1) '.txt'];
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
C = textscan(fid,'%s %*s %f %f %f %*s %f %f %f %*s %f %f %f');
fclose(fid);
Loc = C{1};
Xref = C{2};
Yref = C{3};
Zref = C{4};
Xtar = C{5};
Ytar = C{6};
Ztar = C{7};
Vtx = C{8};
Vty = C{9};
Vtz = C{10};
