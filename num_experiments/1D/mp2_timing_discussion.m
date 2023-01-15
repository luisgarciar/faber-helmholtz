% run mp2_1dsom.m first.

time_fs1 = time_lu_FS(1,:).';
time_fs2 = time_lu_FS(2,:).';
time_fs3 = time_lu_FS(3,:).';

speedup1 = abs(time_fs1 - time_lu)./abs(time_lu);
speedup2 = abs(time_fs2 - time_lu)./abs(time_lu);
speedup3 = abs(time_fs3 - time_lu)./abs(time_lu);

[time_fs1./time_lu, time_fs2./time_lu, time_fs3./time_lu]
