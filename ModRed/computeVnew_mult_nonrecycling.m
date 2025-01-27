%================================================================================
% MOR by using our PMOR algorithm in Proc. ICIAM07
% JJ can be changed according to the accuracy of the reducded model
% is required

%=================================================================================

% the equation is 
% (E0+Rho*cE1) dx/dt +(K0+Kappa*K1+h*K2)x=b U^2(t)/R(T)
% y=C*x; C=[0,...,0,1,0,...0; 0,...,0,1,0,...,0 ]. 
% 1 on the first row is the 12599th element in C, 1 is the 42440th element in C. 


load microhotplate 

% Kapil
n=size(B,1);
x01 = zeros(n,1);
x0_tilde1 = zeros(n,1);
x02 = zeros(n,1);
x0_tilde2 = zeros(n,1);    
num_systems = 0;
U1 = [];
U1_tilde = [];

 JJ=2;  % highest order of the moments included
 
%====================================================================00
% expansion around s0, kappa0, cp0, Rho0, h0

 
 s0=0; kappa0=2.5; cp0=439; Rho0=3100; h0=10;
 
 
 coeff=A0+s0*Rho0*cp0*E1+kappa0*A1+h0*A2;

 
fvnew

V1=V;
  
%====================================================================00
% expansion around s0, kappa0, cp0, Rho0, h0

 
 s0=0; kappa0=3; cp0=500; Rho0=3100; h0=10.5;
 
 
 coeff=A0+s0*Rho0*cp0*E1+kappa0*A1+h0*A2;

 
fvnew

V2=V;
 
%======================================================================
%====================================================================00
% expansion around s0, kappa0, cp0, Rho0, h0

 
 s0=0; kappa0=4; cp0=700; Rho0=3100; h0=11;
 
 
 coeff=A0+s0*Rho0*cp0*E1+kappa0*A1+h0*A2;

 
fvnew

V3=V;
 
%======================================================================

VV=[V1,V2,V3];

V=forthognalize(VV,0.0000001);


n=size(B,1);

Br=V'*B;

e1=E0*V; E0r=V'*e1; clear e1

e1=E1*V; E1r=V'*e1; clear e1

k1=A0*V; A0r=V'*k1; clear k1

k1=A1*V; A1r=V'*k1; clear k1

k1=A2*V; A2r=V'*k1; clear k1

Cr=C*V;

%=============================================================================
%%
figure (1);
hold all;
title('Time Comparison of BiCGSTAB with Recycling BiCGSTAB', ...
    'fontsize', 12, 'fontweight', 'bold');
xlabel('Linear System Number', 'FontSize', 12,'FontWeight', 'bold');
ylabel('Time in Seconds', 'FontSize', 12, 'FontWeight', 'bold');

lengthOfVec = length(vec_bicgstab_time);
plot(1:lengthOfVec, vec_bicgstab_time, '-+r', ...
	1:lengthOfVec, vec_recycling_time, '--ob','Linewidth', 2);

h_legend = legend( ...
    sprintf('BiCGSTAB (Total time %f)', sum(vec_bicgstab_time)), ...
    sprintf('Recycling BiCGSTAB (Total Time %f)', sum(vec_recycling_time)));
set(h_legend, 'FontSize', 12, 'FontWeight', 'bold'); 

%%
figure (2);
hold all;
title('Iteration Comparison of BiCGSTAB with Recycling BiCGSTAB', ...
    'fontsize', 12, 'fontweight', 'bold');
xlabel('Linear System Number', 'FontSize', 12,'FontWeight', 'bold');
ylabel('Number of Iterations', 'FontSize', 12, 'FontWeight', 'bold');

plot(1:lengthOfVec, vec_bicgstab_iter, '-+r', ...
	1:lengthOfVec, vec_recycling_iter, '--ob','Linewidth', 2);

h_legend = legend( ...
    sprintf('BiCGSTAB (Total Iterations %f)', sum(vec_bicgstab_iter)), ...
    sprintf('Recycling BiCGSTAB (Total Iterations %f)', sum(vec_recycling_iter)));
set(h_legend, 'FontSize', 12, 'FontWeight', 'bold'); 

%%
figure (3);
hold all;
title('Residual Comparison of BiCGSTAB with Recycling BiCGSTAB', ...
    'fontsize', 12, 'fontweight', 'bold');
xlabel('Linear System Number', 'FontSize', 12,'FontWeight', 'bold');
ylabel('Residual', 'FontSize', 12, 'FontWeight', 'bold');

plot(1:lengthOfVec, vec_bicgstab_res, '-+r', ...
	1:lengthOfVec, vec_recycling_res, '--ob','Linewidth', 2);

h_legend = legend('BiCGSTAB', 'Recycling BiCGSTAB');
set(h_legend, 'FontSize', 12, 'FontWeight', 'bold'); 




