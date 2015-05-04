function [A_volt_2V,A_current_2V,A_volt_24V,A_current_24V,A_volt_25V,A_current_25V]= Experimental_Data
% Loads the experimental data from http://jes.ecsdl.org/content/152/5/D79

N1_2V = 23;
N2_2V = 51;

N1_24V = 29;
N2_24V = 63;

N1_25V = 31;
N2_25V = 73;

% Experimental Data- 12.7s CC charge
Data_2V = importdata('data_2V.txt');
A2_2V = Data_2V.data;
A3_2V = A2_2V(:,2:3);
A3_2V(1,1) = 0;
A_volt_2V = A3_2V(1:N1_2V,:);
A_current_2V = A3_2V(N1_2V+1:N2_2V,:);

% Experimental Data- 18s CC charge
A = importdata('data.txt');
A2 = A.data;
A3 = A2(:,2:3);
A3(1,1) = 0;
A_volt_24V = A3(1:N1_24V,:);
A_current_24V = A3(N1_24V+1:N2_24V,:);

% Experimental Data- 18.7s CC charge
A_25V = importdata('data_25V.txt');
A2_25V = A_25V.data;
A3_25V = A2_25V(:,2:3);
A3_25V(1,1) = 0;
A_volt_25V = A3_25V(1:N1_25V,:);
A_current_25V = A3_25V(N1_25V+1:N2_25V,:);

end