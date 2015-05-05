% This function simulates a constant current charging profile for the
% supercapacitor model of Drummond, Ross, David A. Howey, and Stephen R. Duncan. "Low-order mathematical modelling of electric double layer supercapacitors using spectral methods." Journal of Power Sources 277 (2015): 317-328.
% The model equations are discretised using a spectral collocation method
% which is implemented using the "cheb.m" function

% The input to the function is the charging time tf. The model outputs are
% the voltages (y_store), currents (i_store), simulation times (t_store),
% model states (states) and otutput system matrices (C_left,D_left,C_right,D_right).
% The model parameters were obtained from the super_params.m file.

% My name is Ross Drummond (ross.drummond@eng.ox.ac.uk) and I hold the MIT license for this code. 
% The accompanying paper for the code can be found at http://www.sciencedirect.com/science/article/pii/S0378775314019739.
% I would ask that you cite this paper as Drummond, Ross, David A. Howey, and Stephen R. Duncan. "Low-order mathematical modelling of electric double layer supercapacitors using spectral methods." Journal of Power Sources 277 (2015): 317-328 if you want to use this code for your own research. 
% For further details on the work of the Energy Power Group at Oxford, please see epg.eng.ox.ac.uk.

function [y_store,i_store,t_store,states,C_left,D_left,C_right,D_right] = Constant_Current(tf)
[Da,Db,Dc,La,Lb,Lc,K1,K2,Kapa_solid,Kapa_elyte,sigma,epsilon_solid,epsilon_elyte,a,C,F,Na,Nb,Nc] = super_params; % Model parameters
[x0,i] = intitial_cons_current(Na,Nb,Nc);% Load the initial conditions and current

%% Spectral collocation Chebyshev differentiation matrices.
[D_cheba,~] = cheb(Na);
[D_chebb,~] = cheb(Nb);
[D_chebc,~] = cheb(Nc);

D_cheba= -D_cheba;
D_chebb= -D_chebb;
D_chebc= -D_chebc;

%% Extract the differentiation matrices for the interior and boundary domains
D_cheba_in = D_cheba(2:Na,2:Na);
D_chebb_in = D_chebb(2:Nb,2:Nb);
D_chebc_in = D_chebc(2:Nc,2:Nc);

DD_cheba = D_cheba^2;
DD_chebb = D_chebb^2;
DD_chebc = D_chebc^2;

DD_cheba_in = DD_cheba(2:Na,2:Na);
DD_chebb_in = DD_chebb(2:Nb,2:Nb);
DD_chebc_in = DD_chebc(2:Nc,2:Nc);

%% Apply the boundary conditions
% BC's concentration
Dc_a = [(2/La)*D_cheba(1,1),(2/La)*D_cheba(1,Na+1), 0, 0;
    (2*Da/La)*D_cheba(Na+1,1), -(2*Db/Lb)*D_chebb(1,1)+(2*Da/La)*D_cheba(Na+1,Na+1), -(2*Db/Lb)*D_chebb(1,Nb+1), 0;
    0, -(2*Db/Lb)*D_chebb(Nb+1,1), (2*Dc/Lc)*D_chebc(1,1)-(2*Db/Lb)*D_chebb(Nb+1,Nb+1),(2*Dc/Lc)*D_chebc(1,Nc+1);
    0, 0, (2/Lc)*D_chebc(Nc+1,1),(2/Lc)*D_chebc(Nc+1,Nc+1)];

Dc_b = [-(2/La)*D_cheba(1,2:Na),zeros(1,Nb-1), zeros(1,Nc-1);
    -(2*Da/La)*D_cheba(Na+1,2:Na),(2*Db/Lb)*D_chebb(1,2:Nb),  zeros(1,Nc-1);
    zeros(1,Na-1), (2*Db/Lb)*D_chebb(Nb+1,2:Nb), -(2*Dc/Lc)*D_chebc(1,2:Nc);
    zeros(1,Na-1), zeros(1,Nb-1), -(2/Lc)*D_chebc(Nc+1,2:Nc)];

%% BC's Phi1
D_volt1_a = [(2/La)*D_cheba(1,1),(2/La)*D_cheba(1,Na+1), 0,0;
    (2/La)*D_cheba(Na+1,1), (2/La)*D_cheba(Na+1,Na+1),0,0;
    0,0, (2/Lc)*D_chebc(1,1),(2/Lc)*D_chebc(1,Nc+1);
    0,0,(2/Lc)*D_chebc(Nc+1,1),(2/Lc)*D_chebc(Nc+1,Nc+1)];

D_volt1_b = [-(2/La)*D_cheba(1,2:Na), zeros(1,Nb-1),zeros(1,Nc-1);
    -(2/La)*D_cheba(Na+1,2:Na), zeros(1,Nb-1),zeros(1,Nc-1);
    zeros(1,Na-1),zeros(1,Nb-1), -(2/Lc)*D_chebc(1,2:Nc);
    zeros(1,Na-1), zeros(1,Nb-1),-(2/Lc)*D_chebc(Nc+1,2:Nc)];

D_volt1_c = [-1/sigma;0; 0 ; -1/sigma];

D_volt1_A = D_volt1_a\D_volt1_b;
D_volt1_B = D_volt1_a\D_volt1_c;

D_volt1_A(1,:) = 0;% Reference: Left hand side Phi1 = 0
D_volt1_B(1,:) = 0;

%% BC's Phi2
D_volt2_a = [(2/La)*D_cheba(1,1), (2/La)*D_cheba(1,Na+1),0,0;
    Kapa_solid*(2/La)*D_cheba(Na+1,1),Kapa_solid*(2/La)*D_cheba(Na+1,Na+1)-Kapa_elyte*(2/Lb)*D_chebb(1,1),-Kapa_elyte*(2/Lb)*D_chebb(1,Nb+1),0;
    0,Kapa_elyte*(2/Lb)*D_chebb(Nb+1,1),-Kapa_solid*(2/Lc)*D_chebc(1,1)+Kapa_elyte*(2/Lb)*D_chebb(Nb+1,Nb+1),-Kapa_solid*(2/Lc)*D_chebc(1,Nc+1);
    0,0,(2/Lc)*D_chebc(Nc+1,1),(2/Lc)*D_chebc(Nc+1,Nc+1);];

D_volt2_b = K2*[0,0,0,0;
    Kapa_solid*(2/La)*D_cheba(Na+1,1),Kapa_solid*(2/La)*D_cheba(Na+1,Na+1)-Kapa_elyte*(2/Lb)*D_chebb(1,1),-Kapa_elyte*(2/Lb)*D_chebb(1,Nb+1),0;
    0,Kapa_elyte*(2/Lb)*D_chebb(Nb+1,1),-Kapa_solid*(2/Lc)*D_chebc(1,1)+Kapa_elyte*(2/Lb)*D_chebb(Nb+1,Nb+1),-Kapa_solid*(2/Lc)*D_chebc(1,Nc+1);
    0,0,0,0;];

D_volt2_c = [-(2/La)*D_cheba(1,2:Na), zeros(1,Nb-1),zeros(1,Nc-1);
    -Kapa_solid*(2/La)*D_cheba(Na+1,2:Na), Kapa_elyte*(2/Lb)*D_chebb(1,2:Nb),zeros(1,Nc-1);
    zeros(1,Na-1), -Kapa_elyte*(2/Lb)*D_chebb(Nb+1,2:Nb),Kapa_solid*(2/Lc)*D_chebc(1,2:Nc);
    zeros(1,Na-1), zeros(1,Nb-1),-(2/Lc)*D_chebc(Nc+1,2:Nc)];

D_volt2_d = K2*[zeros(1,Na-1), zeros(1,Nb-1),zeros(1,Nc-1);
    -Kapa_solid*(2/La)*D_cheba(Na+1,2:Na), Kapa_elyte*(2/Lb)*D_chebb(1,2:Nb),zeros(1,Nc-1);
    zeros(1,Na-1), -Kapa_elyte*(2/Lb)*D_chebb(Nb+1,2:Nb),Kapa_solid*(2/Lc)*D_chebc(1,2:Nc);
    zeros(1,Na-1), zeros(1,Nb-1),zeros(1,Nc-1);];

D_volt2_A= D_volt2_a\(D_volt2_b);
D_volt2_B = D_volt2_a\D_volt2_c;
D_volt2_C = D_volt2_a\D_volt2_d;

% D_volt2_A(1,:) = 0;% Reference: Left hand side Phi2 = 0
% D_volt2_B(1,:) = 0;
% D_volt2_C(1,:) = 0;

%% EQUATION 1- Diffusion Equation
DD_in_eqn1 = zeros(Na+Nb+Nc-3,Na+Nb+Nc-3);
DD_in_eqn1(1:Na-1,1:Na-1) = (4*Da/La^2)*DD_cheba_in;
DD_in_eqn1(Na:Na+Nb-2,Na:Na+Nb-2) = (4*Db/Lb^2)*DD_chebb_in;
DD_in_eqn1(Na+Nb-1:Na+Nb+Nc-3,Na+Nb-1:Na+Nb+Nc-3) = (4*Dc/Lc^2)*DD_chebc_in;

DD_ex_eqn1 = zeros(Na-1+Nb-1+Nc-1,4);
DD_ex_eqn1(1:Na-1,1) = (4*Da/La^2)*DD_cheba(2:Na,1);
DD_ex_eqn1(1:Na-1,2) = (4*Da/La^2)*DD_cheba(2:Na,Na+1);
DD_ex_eqn1(Na:Na+Nb-2,2) = (4*Db/Lb^2)*DD_chebb(2:Nb,1);
DD_ex_eqn1(Na:Na+Nb-2,3) = (4*Db/Lb^2)*DD_chebb(2:Nb,Nb+1);
DD_ex_eqn1(Na+Nb-1:Na+Nb+Nc-3,3) = (4*Dc/Lc^2)*DD_chebc(2:Nc,1);
DD_ex_eqn1(Na+Nb-1:Na+Nb+Nc-3,4) = (4*Dc/Lc^2)*DD_chebc(2:Nc,Nc+1);

D_star = Dc_a\Dc_b;
DD_ex_star = DD_ex_eqn1*D_star;

A_1 = [DD_in_eqn1+DD_ex_star,zeros(Na-1+Nb-1+Nc-1),zeros(Na-1+Nb-1+Nc-1)];
B_1 = zeros(Na-1+Nb-1+Nc-1,1);

%% EQUATION 2- Double Layer Equation
DD_ex_eq2 = zeros(Na-1+Nb-1+Nc-1,4);
DD_ex_eq2(1:Na-1,1) = (4/La^2)*DD_cheba(2:Na,1);
DD_ex_eq2(1:Na-1,2) = (4/La^2)*DD_cheba(2:Na,Na+1);
DD_ex_eq2(Na+Nb-1:Na+Nb+Nc-3,3) = (4/Lc^2)*DD_chebc(2:Nc,1);
DD_ex_eq2(Na+Nb-1:Na+Nb+Nc-3,4) = (4/Lc^2)*DD_chebc(2:Nc,Nc+1);

DD_in_eq2 = zeros(Na+Nb+Nc-3,Na+Nb+Nc-3);
DD_in_eq2(1:Na-1,1:Na-1) = (4/La^2)*DD_cheba_in;
DD_in_eq2(Na:Na+Nb-2,Na:Na+Nb-2) = zeros(Nb-1,Nb-1);
DD_in_eq2(Na+Nb-1:Na+Nb+Nc-3,Na+Nb-1:Na+Nb+Nc-3) = (4/Lc^2)*DD_chebc_in;

A_in_2 = [zeros(Na-1+Nb-1+Nc-1),sigma*DD_in_eq2,zeros(Na-1+Nb-1+Nc-1)];
B_in_2 = zeros(Na-1+Nb-1+Nc-1,1);

A_ex_2 = [zeros(Na-1+Nb-1+Nc-1),sigma*DD_ex_eq2*D_volt1_A ,zeros(Na-1+Nb-1+Nc-1)];
B_ex_2 = sigma*DD_ex_eq2*D_volt1_B;

A_2 = A_in_2+A_ex_2;
B_2 = B_in_2+B_ex_2;

%% EQUATION 3- Algebraic Constaint for curret conservation
D_ex_eq3_volt1 = zeros(Na-1+Nb-1+Nc-1,4);
D_ex_eq3_volt1(1:Na-1,1) = (2/La)*D_cheba(2:Na,1);
D_ex_eq3_volt1(1:Na-1,2) = (2/La)*D_cheba(2:Na,Na+1);
D_ex_eq3_volt1(Na+Nb-1:Na+Nb+Nc-3,3) = (2/Lc)*D_chebc(2:Nc,1);
D_ex_eq3_volt1(Na+Nb-1:Na+Nb+Nc-3,4) = (2/Lc)*D_chebc(2:Nc,Nc+1);
D_ex_eq3_volt1 = sigma*D_ex_eq3_volt1;

D_in_eq3_volt1 = zeros(Na+Nb+Nc-3,Na+Nb+Nc-3);
D_in_eq3_volt1(1:Na-1,1:Na-1) = (2/La)*D_cheba_in;
D_in_eq3_volt1(Na:Na+Nb-2,Na:Na+Nb-2) = zeros(Nb-1,Nb-1);
D_in_eq3_volt1(Na+Nb-1:Na+Nb+Nc-3,Na+Nb-1:Na+Nb+Nc-3) = (2/Lc)*D_chebc_in;
D_in_eq3_volt1= sigma*D_in_eq3_volt1;

D_ex_eq3_volt2 = [Kapa_solid*(2/La)*D_cheba(2:Na,1), Kapa_solid*(2/La)*D_cheba(2:Na,Na+1), zeros(Na-1,1), zeros(Na-1,1);
    zeros(Nb-1,1), Kapa_elyte*(2/Lb)*D_chebb(2:Nb,1), Kapa_elyte*(2/Lb)*D_chebb(2:Nb,Nb+1), zeros(Nb-1,1);
    zeros(Nc-1,1), zeros(Nc-1,1),Kapa_solid*(2/Lc)*D_chebc(2:Nc,1), Kapa_solid*(2/Lc)*D_chebc(2:Nc,Nc+1)];

D_in_eq3_volt2 = zeros(Na+Nb+Nc-3,Na+Nb+Nc-3);
D_in_eq3_volt2(1:Na-1,1:Na-1) = Kapa_solid*(2/La)*D_cheba_in;
D_in_eq3_volt2(Na:Na+Nb-2,Na:Na+Nb-2) = Kapa_elyte*(2/Lb)*D_chebb_in;
D_in_eq3_volt2(Na+Nb-1:Na+Nb+Nc-3,Na+Nb-1:Na+Nb+Nc-3) = Kapa_solid*(2/Lc)*D_chebc_in;

D_ex_eq3_lnc = K2*D_ex_eq3_volt2;
D_in_eq3_lnc= K2*D_in_eq3_volt2;

A_in_3 = [zeros(Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1),D_in_eq3_volt1,D_in_eq3_volt2];
A_in_lnc = D_in_eq3_lnc;
B_in_3 = [ones(Na-1,1);ones(Nb-1,1);ones(Nc-1,1)];

A_ex_3 = [zeros(Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1),D_ex_eq3_volt1*D_volt1_A ,D_ex_eq3_volt2*(D_volt2_B)];
A_ex_lnc = D_ex_eq3_volt2*(D_volt2_C);
A_ex_lnc_ex_BC = -D_ex_eq3_volt2*(D_volt2_A);
A_ex_lnc_ex = D_ex_eq3_lnc;
B_ex_3 = D_ex_eq3_volt1*D_volt1_B;

A_3= A_in_3+A_ex_3;
A_lnc_3 = A_in_lnc +A_ex_lnc;
A_lnextc_3 = A_ex_lnc_ex_BC+A_ex_lnc_ex;

B_3= B_in_3+B_ex_3;

%% A MATRIX
A = zeros(3*(Na-1+Nb-1+Nc-1),3*(Na-1+Nb-1+Nc-1));
A(1:Na-1+Nb-1+Nc-1,:) = A_1;
A(Na-1+Nb-1+Nc:2*(Na-1+Nb-1+Nc-1),:) = A_2;
A(2*(Na-1+Nb-1+Nc-1)+1:3*(Na-1+Nb-1+Nc-1),:) = A_3;

A2 = zeros(3*(Na-1+Nb-1+Nc-1)-Nb+1,3*(Na-1+Nb-1+Nc-1)-Nb+1);
A2(1:Na-1+Nb-1+Nc-1+Na-1,1:Na-1+Nb-1+Nc-1+Na-1) = A(1:Na-1+Nb-1+Nc-1+Na-1,1:Na-1+Nb-1+Nc-1+Na-1);
A2(Na-1+Nb-1+Nc-1+Na:3*(Na-1+Nb-1+Nc-1)-Nb+1,1:Na-1+Nb-1+Nc-1+Na-1) = A(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),1:Na-1+Nb-1+Nc-1+Na-1);

A2(1:Na-1+Nb-1+Nc-1+Na-1,Na-1+Nb-1+Nc-1+Na:3*(Na-1+Nb-1+Nc-1)-Nb+1) = A(1:Na-1+Nb-1+Nc-1+Na-1,Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1));
A2(Na-1+Nb-1+Nc-1+Na:3*(Na-1+Nb-1+Nc-1)-Nb+1,Na-1+Nb-1+Nc-1+Na:3*(Na-1+Nb-1+Nc-1)-Nb+1) = A(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1));

A = A2;
%% A LOG MATRIX
A_lnc = zeros(3*(Na-1+Nb-1+Nc-1),Na-1+Nb-1+Nc-1);
A_lnextc = zeros(3*(Na-1+Nb-1+Nc-1),4);

A_lnc(2*(Na-1+Nb-1+Nc-1)+1:3*(Na-1+Nb-1+Nc-1),1:Na-1+Nb-1+Nc-1) = A_lnc_3;
A_lnextc(2*(Na-1+Nb-1+Nc-1)+1:3*(Na-1+Nb-1+Nc-1),1:4) = A_lnextc_3;

A_lnc= [A_lnc(1:Na-1+Nb-1+Nc-1+Na-1,:);A_lnc(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),:)];
A_lnextc= [A_lnextc(1:Na-1+Nb-1+Nc-1+Na-1,:);A_lnextc(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),:)];

%% B MATRIX
B_2 = [B_2(1:Na-1);B_2(Na-1+Nb:Na-1+Nb-1+Nc-1)];
B = [B_1;B_2;B_3];

%% MASS MATRIX
M = zeros(3*(Na-1+Nb-1+Nc-1),3*(Na-1+Nb-1+Nc-1));
M(1:Na-1+Nb-1+Nc-1,1:Na-1+Nb-1+Nc-1) = epsilon_solid*eye(Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1);% ???????????
M(Na:Na-1+Nb-1,Na:Na-1+Nb-1) = epsilon_elyte*eye(Nb-1,Nb-1);

M(1:Na-1,Na-1+Nb-1+Nc:Na-1+Nb-1+Nc-1+Na-1) = K1*a*C/F*eye(Na-1,Na-1);
M(Na-1+Nb:Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1+Na-1+Nb-1+1:2*(Na-1+Nb-1+Nc-1)) = K1*a*C/F*eye(Nc-1,Nc-1);
M(1:Na-1,2*(Na-1+Nb-1+Nc-1)+1:2*(Na-1+Nb-1+Nc-1)+Na-1) = -K1*a*C/F*eye(Na-1,Na-1);
M(Na-1+Nb:Na-1+Nb-1+Nc-1,2*(Na-1+Nb-1+Nc-1)+Na-1+Nb-1+1:3*(Na-1+Nb-1+Nc-1)) = -K1*a*C/F*eye(Nc-1,Nc-1);

M(Na-1+Nb-1+Nc:2*(Na-1+Nb-1+Nc-1),Na-1+Nb-1+Nc:2*(Na-1+Nb-1+Nc-1)) = a*C*eye(Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1);
M(Na-1+Nb-1+Nc:2*(Na-1+Nb-1+Nc-1),2*(Na-1+Nb-1+Nc-1)+1:3*(Na-1+Nb-1+Nc-1)) = -a*C*eye(Na-1+Nb-1+Nc-1,Na-1+Nb-1+Nc-1);

M = [M(1:Na-1+Nb-1+Nc-1+Na-1,1:Na-1+Nb-1+Nc-1+Na-1),M(1:Na-1+Nb-1+Nc-1+Na-1,Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1));
    M(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),1:Na-1+Nb-1+Nc-1+Na-1),M(Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1),Na-1+Nb-1+Nc-1+Na-1+Nb:3*(Na-1+Nb-1+Nc-1))];

%% SIMULATION
rhsfun = @(t,X) (A*X+A_lnc*log(X(1:Na-1+Nb-1+Nc-1,1))+A_lnextc*log((Dc_a\Dc_b)*X(1:Na-1+Nb-1+Nc-1,1))+B*i); %(define the RHS function)
options = odeset('Mass',M);
[t,states] = ode15s(rhsfun,[0 tf],x0,options); %(solve the problem)
size_states = (size(states,1));

%% OUTPUT MATRICES
D_volt1_out = D_volt1_A ;
B_volt1_out = D_volt1_B;

C_int = D_volt1_out(4,:)-D_volt1_out(1,:);
C_int = [C_int(1:Na-1),C_int(Na-1+Nb:Na-1+Nb-1+Nc-1)];

D = B_volt1_out(4,:)-B_volt1_out(1,:);

C = [zeros(1,Na-1+Nb-1+Nc-1),C_int, zeros(1,Na-1+Nb-1+Nc-1)];

C_left = D_volt1_out(1,:);
C_right = D_volt1_out(4,:);
C_left = [zeros(1,Na-1+Nb-1+Nc-1),C_left(1:Na-1),C_left(Na-1+Nb:Na-1+Nb-1+Nc-1), zeros(1,Na-1+Nb-1+Nc-1)];
C_right = [zeros(1,Na-1+Nb-1+Nc-1),C_right(1:Na-1),C_right(Na-1+Nb:Na-1+Nb-1+Nc-1), zeros(1,Na-1+Nb-1+Nc-1)];

D_left = B_volt1_out(1,:);
D_right = B_volt1_out(4,:);

%% OUTPUT SIGNALS
for j = 1:size_states
    y = -C*states(j,:)'-D*i;
    y_store(j) = y;
    i_store(j) = i;
    t_store(j) = t(j);
end

end

