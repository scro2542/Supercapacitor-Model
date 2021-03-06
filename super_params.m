% The outputs of this function are the model parameters which are defined
% within the code. The parameters were obtained from
% http://jes.ecsdl.org/content/152/5/D79.short.

% My name is Ross Drummond (ross.drummond@eng.ox.ac.uk) and I hold the MIT license for this code. 
% The accompanying paper for the code can be found at http://www.sciencedirect.com/science/article/pii/S0378775314019739.
% I would ask that you cite this paper as Drummond, Ross, David A. Howey, and Stephen R. Duncan. "Low-order mathematical modelling of electric double layer supercapacitors using spectral methods." Journal of Power Sources 277 (2015): 317-328 if you want to use this code for your own research. 
% For further details on the work of the Energy Power Group at Oxford, please see epg.eng.ox.ac.uk.

function [Da,Db,Dc,La,Lb,Lc,K1,K2,Kapa_solid,Kapa_elyte,sigma,epsilon_solid,epsilon_elyte,a,C,F,Na,Nb,Nc] = super_params

Na = 15; % Number of chebyshev nodes in electrode 1
Nb = Na; % Number of chebyshev nodes in the separater
Nc = Na; % Number of chebyshev nodes in electrode 2

epsilon_solid= 0.67; %Porosity constant in the electrode
epsilon_elyte = 0.6; %Porosity constant in the separator

Tortuisity_solid = 2.3; %Tortuosity constant in the electrode
Tortuisity_elyte = 1.29; %Tortuosity constant in the separater

dqp_by_dq = -0.5; %Explained in http://jes.ecsdl.org/content/152/5/D79.
dqn_by_dq = -0.5;

t_neg = 0.45; % Negative transference number
t_pos = 1-t_neg; % Positive transference number

T = 298; %Temperature (K)
F = 9.64853399*10^4; %Faraday costant
c_0 = 930;% Initial concentration (mol/m^3).
C = 1;
a = 42*10^6; %a*C is the specific capacitance (F/m^3)
sigma = 0.0521; %Electrode conductivity

R = 8.314;% Gas constant
f = F/(R*T);% A constant

Kapa_inf = 0.67*10^-1; % Conductivity of the electrolyte in free solution
Kapa_solid= Kapa_inf*epsilon_solid/Tortuisity_solid; % Conductivity of the electrolyte in the electrode
Kapa_elyte= Kapa_inf*epsilon_elyte/Tortuisity_elyte; % Conductivity of the electrolyte in the separater

K1 = (t_neg*dqp_by_dq + t_pos*dqn_by_dq); % A constant
K2 = (t_pos-t_neg)/f; % Another constant

D_solid = 2*Kapa_solid*(t_pos*t_neg/(t_neg+t_pos))*(R*T/((F^2)*c_0));% Diffusion coefficient in the electrode
D_elyte = 2*Kapa_elyte*(t_pos*t_neg/(t_neg+t_pos))*(R*T/((F^2)*c_0));% Diffusion coefficient in the electrolyte

La = 50*10^-6;% Electrode 1 length (m)
Lb = 25*10^-6;% Separater length (m)
Lc = 50*10^-6;% Electrode 2 length (m)

Da =D_solid; % Electrode 1 diffusion co-ef
Db = D_elyte; % Separator diffusion co-ef
Dc = D_solid; % Electrode 1 diffusion co-ef

end
