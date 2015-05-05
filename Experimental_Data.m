% This function loads the experimental data of
% http://jes.ecsdl.org/content/152/5/D79.short against which the model
% simulations are compared. 

% The function outputs are the current and voltage data points for three
% different CC-CV charging profiles. These profiles are differentiated by
% their maximum voltages of 2V, 2.4V and 2.5V.

% My name is Ross Drummond (ross.drummond@eng.ox.ac.uk) and I hold the MIT licencse for this code. 
% The accompanying paper for the code can be found at http://www.sciencedirect.com/science/article/pii/S0378775314019739.
% I would ask that you cite this paper as Drummond, Ross, David A. Howey, and Stephen R. Duncan. "Low-order mathematical modelling of electric double layer supercapacitors using spectral methods." Journal of Power Sources 277 (2015): 317-328 if you want to use this code for your own research. 
% For further details on the work of the Energy Power Group at Oxford, please see epg.eng.ox.ac.uk.

%% The MIT License (MIT)

%Copyright (c) 2015 Ross Drummond

%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.

%%

function [A_volt_2V,A_current_2V,A_volt_24V,A_current_24V,A_volt_25V,A_current_25V]= Experimental_Data

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