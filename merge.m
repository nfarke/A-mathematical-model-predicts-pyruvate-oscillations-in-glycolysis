a1 = load('Resultsup4');
a2 = load('Resultsup5');
a3 = load('Resultsup6');
a4 = load('Resultsup7');
a5 = load('Resultsup8');

IMAG = vertcat(a1.IMAG',a2.IMAG',a3.IMAG',a4.IMAG',a5.IMAG');
PAR = vertcat(a1.PAR,a2.PAR,a3.PAR,a4.PAR,a5.PAR);
RE = vertcat(a1.RE',a2.RE',a3.RE',a4.RE',a5.RE');
Results = vertcat(a1.Results,a2.Results,a3.Results,a4.Results,a5.Results);
stable = vertcat(a1.stable',a2.stable',a3.stable',a4.stable',a5.stable');

save('Results10k','PAR','Results','stable','RE','IMAG','-v7.3')           
