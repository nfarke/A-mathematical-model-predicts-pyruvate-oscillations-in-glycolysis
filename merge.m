a1 = load('Results_1k_5down');
a2 = load('Results_1k_5down2');
a3 = load('Results_1k_5down3');

IMAG = vertcat(a1.IMAG',a2.IMAG',a3.IMAG');
PAR = vertcat(a1.PAR,a2.PAR,a3.PAR);
RE = vertcat(a1.RE',a2.RE',a3.RE');
Results = vertcat(a1.Results,a2.Results,a3.Results);
stable = vertcat(a1.stable',a2.stable',a3.stable');

save('Results_1k_5down_merged','PAR','Results','stable','RE','IMAG','-v7.3')           
