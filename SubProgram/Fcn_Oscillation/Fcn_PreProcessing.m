function Fcn_PreProcessing
% This function is used to calculate some mean matrix coefficients
%
% -------------------------------------
global CI
DuctIndex=CI.setup.DuctIndex;
DuctIndex=[DuctIndex,DuctIndex(end)];
for ss = 1:CI.TP.numSection-1 
    switch DuctIndex(ss+1)
        case 0
            if CI.setup.DuctIndex(ss)==1                                        % From premixer to annular duct
                [G_l_1,G_l_2,F_l] = Fcn_Matrix_calculation_premixer_2_duct(ss); % Calculate Matrices within the premixer section                
                CI.TPM.B1{1,ss}       = G_l_1;                                  % Calculate left-hand side G for the interface -- it is 3X3 matrices
                CI.TPM.B1{2,ss}       = G_l_1;
                CI.TPM.B1_1{1,ss}     = G_l_2;
                CI.TPM.B1_1{2,ss}     = G_l_2;
                CI.TPM.F_l{1,ss}      = F_l;
                CI.TPM.F_l{2,ss}      = F_l;
                
                Theta=CI.TP.Theta(ss);
                if Theta<=1
                    Theta=1.5;
                else
                end
                [B1,B2,B2_inv] = Fcn_Matrix_calculation_WO_Addition_effect(ss, Theta); % Calculate right-hand side G for the interface -- it is 4X4 matrices
                CI.TPM.B2{1,ss}         = B2;
                CI.TPM.B2{2,ss}         = B2;
                CI.TPM.B2_inv{1,ss}     = B2_inv;
                CI.TPM.B2_inv{2,ss}     = B2_inv;
            else
                switch CI.setup.InterfaceIndex(ss+1)
                    case 0
                        Theta=CI.TP.Theta(ss);
                        [B1,B2,B2_inv] = Fcn_Matrix_calculation_WO_Addition_effect(ss,Theta);
                        CI.TPM.B1{1,ss}         = B1;
                        CI.TPM.B1{2,ss}         = B1;
                        CI.TPM.B2{1,ss}         = B2;
                        CI.TPM.B2{2,ss}         = B2;
                        CI.TPM.B2_inv{1,ss}     = B2_inv;
                        CI.TPM.B2_inv{2,ss}     = B2_inv;
                    case {10,11}
                        Theta=CI.TP.Theta(ss);
                        [B1,B2,B2_inv] = Fcn_Matrix_calculation_WO_Addition_effect(ss,Theta);
                        CI.TPM.B1{1,ss}         = B1;               % first step
                        CI.TPM.B1{2,ss}         = B1;
                        CI.TPM.B2{1,ss}         = B2;
                        CI.TPM.B2{2,ss}         = B2;
                        CI.TPM.B2_inv{1,ss}     = B2_inv;
                        CI.TPM.B2_inv{2,ss}     = B2_inv;
                        [B1,B2,B2_inv]= Fcn_Matrix_calculation_W_HA(ss);
                        CI.TPM.B1{1,ss}         = B1;              % second step
                        CI.TPM.B2{1,ss}         = B2;
                        CI.TPM.B2_inv{1,ss}     = B2_inv;
                    otherwise
                        h = msgbox('The interfaces include cases which cannot be considered in the present code.');
                end
            end
        case 1
            if CI.setup.DuctIndex(ss)==0                                        % From annular duct to premixers
                Theta=CI.TP.Theta(ss);
                if Theta>=1
                    Theta=0.5;
                else
                end
                [B1,B2,B2_inv] = Fcn_Matrix_calculation_WO_Addition_effect(ss,Theta); % Calculate left-hand side G for the interface -- it is 4X4 matrices
                CI.TPM.B1{1,ss}         = B1;
                CI.TPM.B1{2,ss}         = B1;
                
                [F_r_inv,G_r_inv] = Fcn_Matrix_calculation_duct_2_premixer(ss); % Calculate Matrices within the premixer section
                CI.TPM.B2_inv{1,ss}     = G_r_inv;                              % Calculate right-hand side G for the interface -- it is 3X3 matrices
                CI.TPM.B2_inv{2,ss}     = G_r_inv;
                CI.TPM.F_r_inv{1,ss}    = F_r_inv;
                CI.TPM.F_r_inv{2,ss}    = F_r_inv;
                
            else if CI.setup.DuctIndex(ss)==1                                   % From premixers to premixers
                    [F_l,G_l,F_r_inv,G_r_inv] = Fcn_Matrix_calculation_premixer_2_premixer(ss);
                    CI.TPM.B1{1,ss}       = G_l;
                    CI.TPM.B1{2,ss}       = G_l;
                    CI.TPM.F_l{1,ss}      = F_l;
                    CI.TPM.F_l{2,ss}      = F_l;
                    CI.TPM.B2_inv{1,ss}   = G_r_inv;
                    CI.TPM.B2_inv{2,ss}   = G_r_inv;
                    CI.TPM.F_r_inv{1,ss}  = F_r_inv;
                    CI.TPM.F_r_inv{2,ss}  = F_r_inv;
             
                else
                    h = msgbox('The ducts include cases which cannot be considered in the present code.');
                end
            end
                    
        otherwise
            h = msgbox('The ducts include cases which cannot be considered in the present code.');
    end
    
    
end
assignin('base','CI',CI)
%
% -------------------------------------------------------------------------
%
function [B1a,B2a,B2a_inv] = Fcn_Matrix_calculation_W_HA(ss)
global CI
u_l=CI.TP.u_mean(2,ss+1);
u_r=CI.TP.u_mean(1,ss+1);
rho_l=CI.TP.rho_mean(2,ss+1);
rho_r=CI.TP.rho_mean(1,ss+1);
p_l=CI.TP.p_mean(2,ss+1);
p_r=CI.TP.p_mean(1,ss+1);
gamma_l=CI.TP.gamma(2,ss+1);
gamma_r=CI.TP.gamma(1,ss+1);
r=CI.setup.R_m;
% ----------------------------------       
B1a(1,1) =0;
B1a(1,2) =u_l;
B1a(1,3) =rho_l;
B1a(1,4) =0;

B1a(2,1) =1;
B1a(2,2) =u_l^2;
B1a(2,3) =2*rho_l*u_l;
B1a(2,4) =0;

B1a(3,1) =0;
B1a(3,2) =0;
B1a(3,3) =0;
B1a(3,4) =r*rho_l*u_l;

B1a(4,1) =gamma_l*u_l/(gamma_l-1);
B1a(4,2) =0.5*u_l^3;
B1a(4,3) =gamma_l/(gamma_l-1)*p_l+1.5*rho_l*u_l^2;
B1a(4,4) =0;

%
B2a(1,1) =0;
B2a(1,2) =u_r;
B2a(1,3) =rho_r;
B2a(1,4) =0;

B2a(2,1) =1;
B2a(2,2) =u_r^2;
B2a(2,3) =2*rho_r*u_r;
B2a(2,4) =0;

B2a(3,1) =0;
B2a(3,2) =0;
B2a(3,3) =0;
B2a(3,4) =r*rho_r*u_r;

B2a(4,1) =gamma_r*u_r/(gamma_r-1);
B2a(4,2) =0.5*u_r^3;
B2a(4,3) =gamma_r/(gamma_r-1)*p_r+1.5*rho_r*u_r^2;
B2a(4,4) =0;

B2a_inv  =[ -(gamma_r*rho_r*u_r^3 - rho_r*u_r^3 + 2*gamma_r*p_r*u_r)/(2*(- rho_r*u_r^2 + gamma_r*p_r)), (gamma_r*p_r - rho_r*u_r^2 + gamma_r*rho_r*u_r^2)/(- rho_r*u_r^2 + gamma_r*p_r),               0, -(rho_r*u_r*(gamma_r - 1))/(- rho_r*u_r^2 + gamma_r*p_r);
            (3*rho_r*u_r^2 - 2*gamma_r*p_r + gamma_r*rho_r*u_r^2)/(2*(rho_r*u_r^3 - gamma_r*p_r*u_r)),                                 (gamma_r*rho_r)/(- rho_r*u_r^2 + gamma_r*p_r),               0,    (rho_r*(gamma_r - 1))/(rho_r*u_r^3 - gamma_r*p_r*u_r);
            (gamma_r*u_r^2 + u_r^2)/(2*(- rho_r*u_r^2 + gamma_r*p_r)),                                  -(gamma_r*u_r)/(- rho_r*u_r^2 + gamma_r*p_r),               0,              (gamma_r - 1)/(- rho_r*u_r^2 + gamma_r*p_r);
            0,                                                                         0, 1/(r*rho_r*u_r),                                                    0];


%
% -------------------------------------------------------------------------
%
function [B1,B2,B2_inv] = Fcn_Matrix_calculation_WO_Addition_effect(ss,Theta)
global CI
Area_l=CI.setup.Area(ss);
Area_r=CI.setup.Area(ss+1);
u_l=CI.TP.u_mean(1,ss);     % Be careful about wheather parameters before or after the flame are used
u_r=CI.TP.u_mean(2,ss+1);
rho_l=CI.TP.rho_mean(1,ss);
rho_r=CI.TP.rho_mean(2,ss+1);
p_l=CI.TP.p_mean(1,ss);
p_r=CI.TP.p_mean(2,ss+1);
gamma_l=CI.TP.gamma(1,ss);
gamma_r=CI.TP.gamma(2,ss+1);
r=CI.setup.R_m;
if Theta>=1
    % Left coefficient Matrix of the conservation equations across the interface
    B1(1,1)=0;
    B1(1,2)=Area_l*u_l;
    B1(1,3)=Area_l*rho_l;
    B1(1,4)=0;                                                       % Coefficients for mass fluctuations on the left side of the interface

    B1(2,1)=Area_r;
    B1(2,2)=Area_l*u_l^2;
    B1(2,3)=2*Area_l*rho_l*u_l;
    B1(2,4)=0;                                                       % Coefficients for axial momentum fluctuations on the left side of the interface

    B1(3,1)=0;
    B1(3,2)=0;
    B1(3,3)=0;
    B1(3,4)=Area_l*r*rho_l*u_l;                                      % Coefficients for circumferential momentum fluctuations on the left side of the interface

    B1(4,1)=gamma_l*u_l*Area_l/(gamma_l-1);
    B1(4,2)=0.5*u_l^3*Area_l;
    B1(4,3)=gamma_l/(gamma_l-1)*p_l*Area_l+1.5*rho_l*u_l^2*Area_l;
    B1(4,4)=0;                                                       % Coefficients for energy fluctuations on the left side of the interface                          
    
    % Right coefficient Matrix of the conservation equations across the interface
    B2(1,1)=0;
    B2(1,2)=Area_r*u_r;
    B2(1,3)=Area_r*rho_r;
    B2(1,4)=0;                                                       % Coefficients for mass fluctuations on the left side of the interface

    B2(2,1)=Area_r;
    B2(2,2)=Area_r*u_r^2;
    B2(2,3)=2*Area_r*rho_r*u_r;
    B2(2,4)=0;                                                       % Coefficients for axial momentum fluctuations on the left side of the interface 

    B2(3,1)=0;
    B2(3,2)=0;
    B2(3,3)=0;
    B2(3,4)=Area_r*r*rho_r*u_r;                                      % Coefficients for circumferential momentum fluctuations on the left side of the interface

    B2(4,1)=gamma_r*u_r*Area_r/(gamma_r-1);
    B2(4,2)=0.5*u_r^3*Area_r;
    B2(4,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;
    B2(4,4)=0;                                                       % Coefficients for energy fluctuations on the left side of the interface
    
    B2_inv =[ -(gamma_r*rho_r*u_r^3 - rho_r*u_r^3 + 2*gamma_r*p_r*u_r)/(2*Area_r*(- rho_r*u_r^2 + gamma_r*p_r)), (gamma_r*p_r - rho_r*u_r^2 + gamma_r*rho_r*u_r^2)/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r)),                   0, -(rho_r*u_r*(gamma_r - 1))/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r));
              (3*rho_r*u_r^2 - 2*gamma_r*p_r + gamma_r*rho_r*u_r^2)/(2*Area_r*(rho_r*u_r^3 - gamma_r*p_r*u_r)),                             (gamma_r*rho_r)/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r)),                   0,    (rho_r*(gamma_r - 1))/(Area_r*(rho_r*u_r^3 - gamma_r*p_r*u_r));
              (gamma_r*u_r^2 + u_r^2)/(2*Area_r*(- rho_r*u_r^2 + gamma_r*p_r)),                              -(gamma_r*u_r)/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r)),                   0,            (gamma_r - 1)/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r));
              0,                                                                         0, 1/(Area_r*r*rho_r*u_r),                                                       0];

    
    
else if Theta<1
    % Left coefficient Matrix of the conservation equations across the interface
    B1(1,1)=0;
    B1(1,2)=Area_l*u_l;
    B1(1,3)=Area_l*rho_l;
    B1(1,4)=0;                                                       % Coefficients for mass fluctuations on the left side of the interface

    B1(2,1)=1/rho_l^(gamma_l);
    B1(2,2)=-gamma_l*p_l/rho_l^(gamma_l+1);
    B1(2,3)=0;
    B1(2,4)=0;                                                       % Coefficients for entropy fluctuations on the left side of the interface

    B1(3,1)=0;
    B1(3,2)=0;
    B1(3,3)=0;
    B1(3,4)=Area_l*r*rho_l*u_l;                                      % Coefficients for circumferential momentum fluctuations on the left side of the interface

    B1(4,1)=gamma_l*u_l*Area_l/(gamma_l-1);
    B1(4,2)=0.5*u_l^3*Area_l;
    B1(4,3)=gamma_l/(gamma_l-1)*p_l*Area_l+1.5*rho_l*u_l^2*Area_l;
    B1(4,4)=0;                                                       % Coefficients for energy fluctuations on the left side of the interface                          
    
    % Right coefficient Matrix of the conservation equations across the interface
    B2(1,1)=0;
    B2(1,2)=Area_r*u_r;
    B2(1,3)=Area_r*rho_r;
    B2(1,4)=0;                                                       % Coefficients for mass fluctuations on the left side of the interface

    B2(2,1)=1/rho_r^(gamma_r);
    B2(2,2)=-gamma_r*p_r/rho_r^(gamma_r+1);
    B2(2,3)=0;
    B2(2,4)=0;                                                       % Coefficients for entropy fluctuations on the left side of the interface 

    B2(3,1)=0;
    B2(3,2)=0;
    B2(3,3)=0;
    B2(3,4)=Area_r*r*rho_r*u_r;                                      % Coefficients for circumferential momentum fluctuations on the left side of the interface

    B2(4,1)=gamma_r*u_r*Area_r/(gamma_r-1);
    B2(4,2)=0.5*u_r^3*Area_r;
    B2(4,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;
    B2(4,4)=0;                                                       % Coefficients for energy fluctuations on the left side of the interface  
    
    B2_inv =[ -(gamma_r*p_r*(2*gamma_r*p_r - 3*rho_r*u_r^2 + 3*gamma_r*rho_r*u_r^2))/(2*Area_r*rho_r*(p_r*gamma_r^2*u_r - rho_r*gamma_r*u_r^3 - p_r*gamma_r*u_r + rho_r*u_r^3)), (rho_r^gamma_r*(gamma_r*p_r - rho_r*u_r^2 + gamma_r*rho_r*u_r^2))/(- p_r*gamma_r^2 + rho_r*gamma_r*u_r^2 + p_r*gamma_r - rho_r*u_r^2),                   0, -(gamma_r*p_r)/(Area_r*(rho_r*u_r^3 - gamma_r*p_r*u_r));
              -(2*gamma_r*p_r - 3*rho_r*u_r^2 + 3*gamma_r*rho_r*u_r^2)/(2*Area_r*(p_r*gamma_r^2*u_r - rho_r*gamma_r*u_r^3 - p_r*gamma_r*u_r + rho_r*u_r^3)),                               (gamma_r*rho_r*rho_r^gamma_r)/(- p_r*gamma_r^2 + rho_r*gamma_r*u_r^2 + p_r*gamma_r - rho_r*u_r^2),                   0,       -rho_r/(Area_r*(rho_r*u_r^3 - gamma_r*p_r*u_r));
              -(2*p_r*gamma_r^2 + rho_r*gamma_r*u_r^2 - rho_r*u_r^2)/(2*Area_r*rho_r*(- p_r*gamma_r^2 + rho_r*gamma_r*u_r^2 + p_r*gamma_r - rho_r*u_r^2)),                                -(gamma_r*rho_r^gamma_r*u_r)/(- p_r*gamma_r^2 + rho_r*gamma_r*u_r^2 + p_r*gamma_r - rho_r*u_r^2),                   0,           -1/(Area_r*(- rho_r*u_r^2 + gamma_r*p_r));
              0,                                                                                                             0, 1/(Area_r*r*rho_r*u_r),                                             0]; 
    end
end
%
% -------------------------------------------------------------------------
%
function [F_r_inv,G_r_inv] = Fcn_Matrix_calculation_duct_2_premixer(ss)
global CI
Area_r=CI.setup.Area(ss+1);   % "_r" denotes the parameters in the premixer section
u_r=CI.TP.u_mean(2,ss+1);
rho_r=CI.TP.rho_mean(2,ss+1);
p_r=CI.TP.p_mean(2,ss+1);
gamma_r=CI.TP.gamma(2,ss+1);
c_r=CI.TP.c_mean(2,ss+1);
%
G_r(1,1)=0;
G_r(1,2)=Area_r*u_r;
G_r(1,3)=Area_r*rho_r;

G_r(2,1)=1/rho_r^(gamma_r);
G_r(2,2)=-gamma_r*p_r/rho_r^(gamma_r+1);
G_r(2,3)=0;

G_r(3,1)=gamma_r*u_r*Area_r/(gamma_r-1);
G_r(3,2)=0.5*u_r^3*Area_r;
G_r(3,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;

G_r_inv=inv(G_r);

F_r(1,1)=1;
F_r(1,2)=1;
F_r(1,3)=0;

F_r(2,1)=1/c_r^2;
F_r(2,2)=1/c_r^2;
F_r(2,3)=-1/c_r^2;

F_r(3,1)=1/rho_r/c_r;
F_r(3,2)=-1/rho_r/c_r;
F_r(3,3)=0;

F_r_inv=inv(F_r);

%
% -------------------------------------------------------------------------
%
function [F_l,G_l,F_r_inv,G_r_inv] = Fcn_Matrix_calculation_premixer_2_premixer(ss)
global CI
Area_l=CI.setup.Area(ss);   % "_l" denotes the parameters in the left section
u_l=CI.TP.u_mean(2,ss);
rho_l=CI.TP.rho_mean(2,ss);
p_l=CI.TP.p_mean(2,ss);
gamma_l=CI.TP.gamma(2,ss);
c_l=CI.TP.c_mean(2,ss);
%
Area_r=CI.setup.Area(ss+1);   % "_r" denotes the parameters in the right section
u_r=CI.TP.u_mean(2,ss+1);
rho_r=CI.TP.rho_mean(2,ss+1);
p_r=CI.TP.p_mean(2,ss+1);
gamma_r=CI.TP.gamma(2,ss+1);
c_r=CI.TP.c_mean(2,ss+1);

Theta=CI.TP.Theta(ss);      % This is just Area_r/Area_l

if Theta>=1
    G_l(1,1)=0;
    G_l(1,2)=Area_l*u_l;
    G_l(1,3)=Area_l*rho_l;
    G_l(2,1)=Area_r;
    G_l(2,2)=Area_l*u_l^2;
    G_l(2,3)=2*Area_l*rho_l*u_l;
    G_l(3,1)=gamma_l*u_l*Area_l/(gamma_l-1);
    G_l(3,2)=0.5*u_l^3*Area_l;
    G_l(3,3)=gamma_l/(gamma_l-1)*p_l*Area_l+1.5*rho_l*u_l^2*Area_l;

    F_l(1,1)=1;
    F_l(1,2)=1;
    F_l(1,3)=0;
    F_l(2,1)=1/c_l^2;
    F_l(2,2)=1/c_l^2;
    F_l(2,3)=-1/c_l^2;
    F_l(3,1)=1/rho_l/c_l;
    F_l(3,2)=-1/rho_l/c_l;
    F_l(3,3)=0;
    
    G_r(1,1)=0;
    G_r(1,2)=Area_r*u_r;
    G_r(1,3)=Area_r*rho_r;
    G_r(2,1)=Area_r;
    G_r(2,2)=Area_r*u_r^2;
    G_r(2,3)=2*Area_r*rho_r*u_r;
    G_r(3,1)=gamma_r*u_r*Area_r/(gamma_r-1);
    G_r(3,2)=0.5*u_r^3*Area_r;
    G_r(3,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;
    G_r_inv=inv(G_r);

    F_r(1,1)=1;
    F_r(1,2)=1;
    F_r(1,3)=0;
    F_r(2,1)=1/c_r^2;
    F_r(2,2)=1/c_r^2;
    F_r(2,3)=-1/c_r^2;
    F_r(3,1)=1/rho_r/c_r;
    F_r(3,2)=-1/rho_r/c_r;
    F_r(3,3)=0;
    F_r_inv=inv(F_r);  
else
    G_l(1,1)=0;
    G_l(1,2)=Area_l*u_l;
    G_l(1,3)=Area_l*rho_l;
    G_l(2,1)=1/rho_l^(gamma_l);
    G_l(2,2)=-gamma_l*p_l/rho_l^(gamma_l+1);
    G_l(2,3)=0;
    G_l(3,1)=gamma_l*u_l*Area_l/(gamma_l-1);
    G_l(3,2)=0.5*u_l^3*Area_l;
    G_l(3,3)=gamma_l/(gamma_l-1)*p_l*Area_l+1.5*rho_l*u_l^2*Area_l;

    F_l(1,1)=1;
    F_l(1,2)=1;
    F_l(1,3)=0;
    F_l(2,1)=1/c_l^2;
    F_l(2,2)=1/c_l^2;
    F_l(2,3)=-1/c_l^2;
    F_l(3,1)=1/rho_l/c_l;
    F_l(3,2)=-1/rho_l/c_l;
    F_l(3,3)=0;
    
    G_r(1,1)=0;
    G_r(1,2)=Area_r*u_r;
    G_r(1,3)=Area_r*rho_r;
    G_r(2,1)=1/rho_r^(gamma_r);
    G_r(2,2)=-gamma_r*p_r/rho_r^(gamma_r+1);
    G_r(2,3)=0;
    G_r(3,1)=gamma_r*u_r*Area_r/(gamma_r-1);
    G_r(3,2)=0.5*u_r^3*Area_r;
    G_r(3,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;
    G_r_inv=inv(G_r);

    F_r(1,1)=1;
    F_r(1,2)=1;
    F_r(1,3)=0;
    F_r(2,1)=1/c_r^2;
    F_r(2,2)=1/c_r^2;
    F_r(2,3)=-1/c_r^2;
    F_r(3,1)=1/rho_r/c_r;
    F_r(3,2)=-1/rho_r/c_r;
    F_r(3,3)=0;
    F_r_inv=inv(F_r);
end

%
% -------------------------------------------------------------------------
%
function [G_l_1,G_l_2,F_l] = Fcn_Matrix_calculation_premixer_2_duct(ss)
global CI
Area_r=CI.setup.Area(ss);   % "_r" denotes the parameters in the premixer section
u_r=CI.TP.u_mean(2,ss);
rho_r=CI.TP.rho_mean(2,ss);
p_r=CI.TP.p_mean(2,ss);
gamma_r=CI.TP.gamma(2,ss);
c_r=CI.TP.c_mean(2,ss);
%
G_l(1,1)=0;
G_l(1,2)=Area_r*u_r;
G_l(1,3)=Area_r*rho_r;

G_l(2,1)=1/rho_r^(gamma_r);
G_l(2,2)=-gamma_r*p_r/rho_r^(gamma_r+1);
G_l(2,3)=0;

G_l(3,1)=gamma_r*u_r*Area_r/(gamma_r-1);
G_l(3,2)=0.5*u_r^3*Area_r;
G_l(3,3)=gamma_r/(gamma_r-1)*p_r*Area_r+1.5*rho_r*u_r^2*Area_r;


F_l(1,1)=1;
F_l(1,2)=1;
F_l(1,3)=0;

F_l(2,1)=1/c_r^2;
F_l(2,2)=1/c_r^2;
F_l(2,3)=-1/c_r^2;

F_l(3,1)=1/rho_r/c_r;
F_l(3,2)=-1/rho_r/c_r;
F_l(3,3)=0;

% Note that G_r is related to the circumferential wavenumber n and thus
% needs to be further caluclated in Fcn_DetEqn_Linear.m or Fcn_DetEqn_NonLinear.m
G_l_1     =G_l;
G_l_1(2,1)=0;
G_l_1(2,2)=Area_r*u_r^2;
G_l_1(2,3)=2*Area_r*rho_r*u_r;
G_l_2     =zeros(3,3);
G_l_2(2,1)=CI.setup.Area(ss+1);

%
%----------------------------------end-------------------------------------