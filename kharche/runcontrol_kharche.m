T = 10000.0;
[vs,is,ts,peak_times] = kharche_SA(T,[],[],1e-12);
save(['kharche_control.mat'],'vs','is','ts','peak_times');



