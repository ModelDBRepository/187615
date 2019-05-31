T = 20000.0;
ts_control = 0.1:0.1:T;
[vs,is,ts,peak_times] = severi_SA(T,0);
vs = interpolate(ts,vs',ts_control)';
is = interpolate_multidim(ts,is',ts_control)';
ts = ts_control;
save(['severi_control.mat'],'vs','is','ts','peak_times');



