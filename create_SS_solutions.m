%this file generates the .mat-file PAR.mat. It generates steady state
%solutions of three models (base strain model, knock-out model, 2xCra model)
rng(15)
ParSize = 2000; %adjust number of Parameter sets
tspan = 0:500;
par_end = zeros(26,ParSize);
opts = odeset('RelTol',1e-08,'AbsTol',1E-12);

%syms Vmax1 Vmax2 Vmax3 Vmax4 Vmax5 Vmax6 k1 k2 k3 Km1 Km2 Km3 Km4 Km5 n1 n2 n3 ratio1 a1 a2 a3 a4 g6p fbp pep pyr
%parpool(18)
cntr = 1;
parfor i1 = 1:ParSize
     disp(i1)
     value = 0; %repeat when stability criteria are not met
     [p,par] = Sample(1,1); %samples one parameter set from case gg = 1
     y0 = par(end-3:end); %initial conditions
            
        %calculate constraint parameters 


        k1 = par(find(strcmp(p,'k1')),1);
        k2 = par(find(strcmp(p,'k2')),1);
        k3 = par(find(strcmp(p,'k3')),1);
        Km1 = par(find(strcmp(p,'Km1')),1);
        Km2 = par(find(strcmp(p,'Km2')),1);
        Km3 = par(find(strcmp(p,'Km3')),1);
        Km4 = par(find(strcmp(p,'Km4')),1);
        Km5 = par(find(strcmp(p,'Km5')),1);

        n1 = par(find(strcmp(p,'n1')),1);
        n2 = par(find(strcmp(p,'n2')),1);
        n3 = par(find(strcmp(p,'n3')),1);
        a1 = par(find(strcmp(p,'a1')),1);
        a2 = par(find(strcmp(p,'a2')),1);
        a3 = par(find(strcmp(p,'a3')),1);

        ratio1 = par(find(strcmp(p,'ratio1')),1);
        x1 = -(ratio1*66.67/10)/(ratio1 - 1);

        g6p = par(find(strcmp(p,'g6p')),1);
        fbp = par(find(strcmp(p,'fbp')),1);
        pep = par(find(strcmp(p,'pep')),1);
        pyr = par(find(strcmp(p,'pyr')),1);

        par(strcmp('Vmax1',p)) = 66.67/10*(k1*pyr/pep + k2 + k3*pyr/pep + 1);
        par(strcmp('Vmax2',p)) = (66.67/10+x1)*(1 + (Km1/g6p).^n1) * pep^a1;
        par(strcmp('Vmax3',p)) = x1*(1 + (Km2/fbp))/pep^a2;
        par(strcmp('Vmax4',p)) = 66.67/10 * (1 +  Km3/fbp);
        par(strcmp('Vmax5',p)) = 66.67/10 * (1 + (Km4/pep)^n2)/fbp^a3;
        par(strcmp('Vmax6',p)) = 66.67/10 * (1 + (Km5/pyr)^n3);

        [~,y]  =  ode23tb(@(t,c) odemodel(t,c,p,par),tspan,y0,opts);  

        %Variables at t(end)
        f6p  =  y(20,1);
        fbp =  y(20,2);
        pep =  y(20,3);
        pyr =  y(20,4);   

        Vmax1 = par(find(strcmp(p,'Vmax1')),1);
        Vmax2 = par(find(strcmp(p,'Vmax2')),1);
        Vmax3 = par(find(strcmp(p,'Vmax3')),1);
        Vmax4 = par(find(strcmp(p,'Vmax4')),1);
        Vmax5 = par(find(strcmp(p,'Vmax5')),1);
        Vmax6 = par(find(strcmp(p,'Vmax6')),1);

        r1 = Vmax1 /((k1*pyr/pep + k2 + k3*pyr/pep + 1));
        r2 = Vmax2 * 1/(1 + (Km1/g6p)^n1) * pep^-a1;
        r3 = Vmax3 * fbp/(fbp + Km2) * pep^a2;
        r4 = Vmax4 * fbp/(fbp + Km3);
        r5 = Vmax5 * 1/(1 + (Km4/pep)^n2) * fbp^a3;
        r6 = Vmax6 * 1/(1 + (Km5/pyr)^n3);

        
        J = [-(Km1*Vmax2*n1*(Km1/g6p)^(n1 - 1))/(g6p^2*pep^a1*((Km1/g6p)^n1 + 1)^2),                                                 (Vmax3*pep^a2)/(Km2 + fbp) - (Vmax3*fbp*pep^a2)/(Km2 + fbp)^2, (Vmax1*((k1*pyr)/pep^2 + (k3*pyr)/pep^2))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2 + (Vmax2*a1)/(pep^(a1 + 1)*((Km1/g6p)^n1 + 1)) + (Vmax3*a2*fbp*pep^(a2 - 1))/(Km2 + fbp),                                                        -(Vmax1*(k1/pep + k3/pep))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2;
              (Km1*Vmax2*n1*(Km1/g6p)^(n1 - 1))/(g6p^2*pep^a1*((Km1/g6p)^n1 + 1)^2), (Vmax4*fbp)/(Km3 + fbp)^2 - (Vmax3*pep^a2)/(Km2 + fbp) - Vmax4/(Km3 + fbp) + (Vmax3*fbp*pep^a2)/(Km2 + fbp)^2,                                                                                    - (Vmax2*a1)/(pep^(a1 + 1)*((Km1/g6p)^n1 + 1)) - (Vmax3*a2*fbp*pep^(a2 - 1))/(Km2 + fbp),                                                                                                                          0;
                                                                     0,              (2*Vmax4)/(Km3 + fbp) - (2*Vmax4*fbp)/(Km3 + fbp)^2 - (Vmax5*a3*fbp^(a3 - 1))/((Km4/pep)^n2 + 1),                - (Vmax1*((k1*pyr)/pep^2 + (k3*pyr)/pep^2))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2 - (Km4*Vmax5*fbp^a3*n2*(Km4/pep)^(n2 - 1))/(pep^2*((Km4/pep)^n2 + 1)^2),                                                                      (Vmax1*(k1/pep + k3/pep))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2;
                                                                     0,                                                                    (Vmax5*a3*fbp^(a3 - 1))/((Km4/pep)^n2 + 1),                  (Vmax1*((k1*pyr)/pep^2 + (k3*pyr)/pep^2))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2 + (Km4*Vmax5*fbp^a3*n2*(Km4/pep)^(n2 - 1))/(pep^2*((Km4/pep)^n2 + 1)^2), - (Vmax1*(k1/pep + k3/pep))/(k2 + (k1*pyr)/pep + (k3*pyr)/pep + 1)^2 - (2*Km5*Vmax6*n3*(Km5/pyr)^(n3 - 1))/(pyr^2*((Km5/pyr)^n3 + 1)^2)];

        EV = eig(J);
        [realx,idx] = max(real(EV));
        imagx = imag(EV);
        imagx = imagx(idx);
        
        RE(i1) = realx;
        IMAG(i1) = imagx;
        ratio(i1) = abs(imagx)/abs(realx);
        
        if realx > -1E-05
           stable(i1) = 0;
        else
           stable(i1) = 1;
        end
        
        pyr =  y(:,4);   
        Results(i1,:) = pyr;
        PAR(i1,:) = par;         
end
                       
save('Resultsup8','PAR','Results','stable','RE','IMAG','-v7.3')           

function dcdt  =  odemodel(t,c,p,par,~)

par(end-3:end,1) = c;

%get kinetic parameters
Vmax1 = par(find(strcmp(p,'Vmax1')),1);
Vmax2 = par(find(strcmp(p,'Vmax2')),1);
Vmax3 = par(find(strcmp(p,'Vmax3')),1);
Vmax4 = par(find(strcmp(p,'Vmax4')),1);
Vmax5 = par(find(strcmp(p,'Vmax5')),1);
Vmax6 = par(find(strcmp(p,'Vmax6')),1);

if t>50
   Vmax1 = Vmax1 * 10;
end

k1 = par(find(strcmp(p,'k1')),1);
k2 = par(find(strcmp(p,'k2')),1);
k3 = par(find(strcmp(p,'k3')),1);
Km1 = par(find(strcmp(p,'Km1')),1);
Km2 = par(find(strcmp(p,'Km2')),1);
Km3 = par(find(strcmp(p,'Km3')),1);
Km4 = par(find(strcmp(p,'Km4')),1);
Km5 = par(find(strcmp(p,'Km5')),1);

n1 = par(find(strcmp(p,'n1')),1);
n2 = par(find(strcmp(p,'n2')),1);
n3 = par(find(strcmp(p,'n3')),1);
a1 = par(find(strcmp(p,'a1')),1);
a2 = par(find(strcmp(p,'a2')),1);
a3 = par(find(strcmp(p,'a3')),1);

g6p = par(find(strcmp(p,'g6p')),1);
fbp = par(find(strcmp(p,'fbp')),1);
pep = par(find(strcmp(p,'pep')),1);
pyr = par(find(strcmp(p,'pyr')),1);

r1 = Vmax1 /((k1*pyr/pep + k2 + k3*pyr/pep + 1));
r2 = Vmax2 * 1/(1 + (Km1/g6p)^n1) * pep^-a1;
r3 = Vmax3 * fbp/(fbp + Km2) * pep^a2;
r4 = Vmax4 * fbp/(fbp + Km3);
r5 = Vmax5 * 1/(1 + (Km4/pep)^n2) * fbp^a3;
r6 = Vmax6 * 1/(1 + (Km5/pyr)^n3);


dcdt(1,1) = r1 + r3 - r2; %g6p
dcdt(2,1) = r2 - r3 - r4; %fbp
dcdt(3,1) = 2*r4 - r5 - r1;%pep
dcdt(4,1) = r1 + r5 - 2*r6;%pyr

end

