function [x0,i] = intitial_cons_current(Na,Nb,Nc)
% The simulation initial conditions x0 and charging current 
i = 100/2.747;% Constant current magnitude (A)

%% INITIAL CONDITIONS
% Concentraion
c_0 = zeros(Na-1+Nb-1+Nc-1,1);
c_0(1:Na-1)= 0.93*1000;
c_0(Na:Na-1+Nb-1)= 0.93*1000;
c_0(Na-1+Nb:Na-1+Nb-1+Nc-1)= 0.93*1000;

%% Phi1
[~,x_CHEBa] = cheb(Na);
[~,x_CHEBb] = cheb(Nb);
[~,x_CHEBc] = cheb(Nc);

La = 50*10^-6;
Lb = 25*10^-6;
Lc = 50*10^-6;

sigma = 0.0521;

x = zeros(1,Na+1+(Nb-2)+Nc+1);
x(1:Na+1) = La/2*(x_CHEBa(Na+1:-1:1)+1);
x(Na+2:Na+2+Nb-2) = La+Lb/2*(x_CHEBb(Nb:-1:2)+1);
x(Na+2+Nb-1:Na+2+Nb-1+Nc) = La+Lb+Lc/2*(x_CHEBc(Nc+1:-1:1)+1);

volt1_0 = zeros(Na-1+Nc-1,1);
volt1_0(1:Na-1)= -i*x(1:Na-1)/sigma;
volt1_0(Na:Na-1+Nc-1)= -i*x(Na+Nb-1:Na+Nb+Nc-3)/sigma-1.5;

%% Phi2
volt2_0 = zeros(Na-1+Nb-1+Nc-1,1);
v2 = 0.1;

volt2_0(1:Na-1)= v2;
volt2_0(Na:Na-1+Nb-1)= v2;
volt2_0(Na-1+Nb:Na-1+Nb-1+Nc-1)= v2;

%%
x0 = [c_0;volt1_0; volt2_0];

end

