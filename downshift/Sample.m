function [p,par]=Sample(ParSize,gg)

p={'Vmax1','Vmax2','Vmax3','Vmax4','Vmax5','Vmax6', 'k1','k2','k3','Km1','Km2',...
    'Km3','Km4','Km5','n1','n2','n3','ratio1','a1','a2','a3','a4','g6p','fbp','pep','pyr'};

%parameters average flux
if gg == 1 %called by create_SS_solutions_par -> for figure4b
    par_lb=[1,1,1,1,1,1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 0.01, 1,1,1,1   1,1,1,1]';
    par_ub=[1,1,1,1,1,1, 10,  10,  10,  10,  10,  10,  10,  10, 4, 4, 4,  1,     4,4,4,4  1,1,1,1]';
end

par=zeros(length(par_lb),ParSize);
%random sampling
for i=1:ParSize
    par(:,i)= 10.^((log10(par_lb)-log10(par_ub)).*rand(length(par_lb), 1)+log10(par_ub));
end

end