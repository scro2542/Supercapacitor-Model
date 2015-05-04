clc; clear; close all;

% This file is the master file for a constant-current(CC), constant-voltage (CV) charging profile
% for the Verbrugge supercapacitor model discretised using the spectral
% collocation method. The simulation results are compared to experimental data. Both the model parameters and experimental data were obtained from http://jes.ecsdl.org/content/152/5/D79.short.

% In the CC charging profile, the model input is the current and the
% output is the voltage. With the CV charging profile, the input is the
% voltage and the current is the output. These simulations are run in
% Constant_Current.m and Constant_voltage.m. The spectral collocation
% differentiation matrices were obtained using cheb.m from
% https://people.maths.ox.ac.uk/trefethen/spectral.html.  The model parameters are set up in super_params.m and are defined there.

% To simulate the results of http://jes.ecsdl.org/content/152/5/D79.short, the 3 CC charging times (12.7s,18.7s
% and 23.2s) were implemented. The experimental results for these charging profiles are stored in different files data_2V.txt, data.txt and data_25V.txt. 

% My name is Ross Drummond (ross.drummond@eng.ox.ac.uk) and I hold the liscence for this code. The
% accompanying paper for the code can be found at http://www.sciencedirect.com/science/article/pii/S0378775314019739.
% I would ask that you cite this paper if you want to use this code for your own research. For further details on the
% work of the Energy Power Group at Oxford, please see epg.eng.ox.ac.uk.

%% Model Set-up
[Da,Db,Dc,La,Lb,Lc,K1,K2,Kapa_solid,Kapa_elyte,sigma,epsilon_solid,epsilon_elyte,a,C,F,Na,Nb,Nc] = super_params; %The model parameters

[D_cheba,x_CHEBa] = cheb(Na); [D_chebb,x_CHEBb] = cheb(Nb); [D_chebc,x_CHEBc] = cheb(Nc);% Spectral collocation differentiation matrices.

x = zeros(1,Na+1+(Nb-2)+Nc+1); %Converting the Chebyshev spatial distribution to the true domain 
x(1:Na+1) = La/2*(x_CHEBa(Na+1:-1:1)+1);
x(Na+2:Na+2+Nb-2) = La+Lb/2*(x_CHEBb(Nb:-1:2)+1);
x(Na+2+Nb-1:Na+2+Nb-1+Nc) = La+Lb+Lc/2*(x_CHEBc(Nc+1:-1:1)+1);

[~,i] = intitial_cons_current(Na,Nb,Nc);% Load the current

% tf_C = 12.7;% CC Charging time for 2V max volt charge
% tf_C = 18;% CC Charging time for 2.4V max volt charge
tf_C = 23.2;% CC Charging time for 2.5V max volt charge

%% Constant Current (CC) Charge
[y_store_C,i_store_C,t_store_C,states_C,C_left,D_left,C_right,D_right] = Constant_Current(tf_C);% Simulates the CC charging profiles
size_states_C = (size(states_C,1)); 

%% Constant Potential (CV) Charge
tf_V= 6; % Time for the CV charging profile
x0 = states_C(size_states_C,:)'; % Obtain initial conditions for the CV simulation from the end state values of CC charging profile

V_left = C_left*x0+D_left*i; % Obtain the potential at the left current collector 
V_right = C_right*x0+D_right*i; % Obtain the potential at the right current collector

V_right = V_right+0.675; % Perturb these potentials to form a voltage jump
V_right = V_right+0.175;
V_right = V_right+0.22;
i_V = [V_left;V_right]; %New input vector is the voltage

[y_store_V,i_store_V,t_store_V,states_V] = Constant_voltage(tf_V,x0,i_V);% Simulates the CV charging profile
size_states_V = (size(states_V,1));

%% Post Processing
t_store_V_aug = t_store_V+ tf_C;

i_store_V = 2.747*i_store_V; %Convert from current density (A/m^2) to current (A)
i_store_C = 2.747*i_store_C; % Current collector surface area is 2.747m^2

t_store = [t_store_C';t_store_V_aug'];% Combine CC and CV charging simulations
y_store_Volt = [y_store_C';y_store_V'];
y_store_Curr = [i_store_C';i_store_V'];
states_aug = [states_C;states_V];

[A_volt_2V,A_current_2V,A_volt_24V,A_current_24V,A_volt_25V,A_current_25V]= Experimental_Data; % Load experimental data from Verburgge

%% Plotting
font_size = 12;

h_25V = figure;% Plots the simulations results and experimental data for the CC CV charging profile
h3 = subplot(2,1,1);
plot(h3,t_store,y_store_Volt,'-k',A_volt_25V(:,1),A_volt_25V(:,2),'xb','LineWidth',3,'markersize',12)
xlabel('Time (s)','interpreter','latex', 'FontSize', font_size )
ylabel('Voltage (V)','interpreter','latex', 'FontSize', font_size )
legend('Numerical Simulation','Experimental Data','Location','southwest')
grid on
set(gca,'FontSize',8)
h2 = subplot(2,1,2);
plot(h2,t_store,y_store_Curr,'-k',A_current_25V(:,1),A_current_25V(:,2),'xb','LineWidth',3,'markersize',12)
axis([0 max(t_store_V_aug) -1200 200]);
xlabel('Time (s)','interpreter','latex', 'FontSize',font_size )
ylabel('Current (A)','interpreter','latex', 'FontSize', font_size )
legend('Numerical Simulation','Experimental Data','Location','southwest')
grid on
set(gca,'FontSize',8)

%% Plots the evolution of the states along the charging profiles

[X,Y] = meshgrid([x(2:Na),x(Na+2:Na+Nb),x(Na+Nb+2:Na+Nb+Nc)],t_store);

figure % Concentraion
surf(X,Y,states_aug(:,1:Na+Nb+Nc-3),'EdgeColor','none')
xlabel('x(m)','interpreter','latex', 'FontSize', font_size');
ylabel('t(s)','interpreter','latex', 'FontSize', font_size');
zlabel('c (mol/m$^3$)','interpreter','latex', 'FontSize', font_size')
title('Variation of c','interpreter','latex', 'FontSize', font_size');

volt1 = [states_aug(:,Na+Nb+Nc-3+1:Na+Nb+Nc-3+Na-1),zeros(max(size(t_store)),Nb-1),states_aug(:,Na+Nb+Nc-3+Na-1+1:Na+Nb+Nc-3+Na-1+Nc-1)];
figure%Phi 1
surf(X,Y,volt1,'EdgeColor','none')
xlabel('x(m)','interpreter','latex', 'FontSize', font_size');
ylabel('t(s)','interpreter','latex', 'FontSize', font_size');
zlabel('$\phi_1$','interpreter','latex', 'FontSize', font_size')
title('Variation of $\phi_1$','interpreter','latex', 'FontSize', font_size');
colorbar

figure %Phi2
surf(X,Y,states_aug(:,2*(Na+Nb+Nc-3)+1-Nb+1:3*(Na+Nb+Nc-3)-Nb+1),'EdgeColor','none')
xlabel('x (m)','interpreter','latex', 'FontSize', font_size');
ylabel('t (s)','interpreter','latex', 'FontSize', font_size');
zlabel('$\phi_2$','interpreter','latex', 'FontSize', font_size');
title('Variation of $\Phi_2$','interpreter','latex', 'FontSize', font_size');
colorbar

