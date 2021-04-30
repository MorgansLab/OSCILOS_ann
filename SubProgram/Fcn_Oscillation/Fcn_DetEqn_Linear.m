function F = Fcn_DetEqn_Linear(s)
% This function only accounts for linear cases

omega=-1i*s;       % This depends on the defination of s.

global CI
n=CI.setup.n;
r=CI.setup.R_m;
if n==0
    gamma_r_n=1;
else
    gamma_r_n=sin(n*pi/CI.setup.Premixer_Num)/(n*pi/CI.setup.Premixer_Num);
end

%% Boundary condition
% Inlet
M=CI.TP.M_mean(1,1);
c=CI.TP.c_mean(1,1);
u=CI.TP.u_mean(1,1);
gamma=CI.TP.gamma(1,1);
% Defination of k and alpha in the duct
k_p=(M*omega-(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
k_n=(M*omega+(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
k_0=-omega/u;
alpha_p=omega+u*k_p;
alpha_n=omega+u*k_n;
if strcmp(CI.setup.inlet,'open')==1
    R1=-1;
    R_E_1=0;
    R_V_1=0;
else if strcmp(CI.setup.inlet,'closed')==1
        R1=-k_n/alpha_n/(k_p/alpha_p);
        R_E_1=0;
        R_V_1=0;
    else if strcmp(CI.setup.inlet,'choked')==1 % needs validation
            R1   =-(c*k_n/alpha_n+n^2*c/r^2/k_0/alpha_n-gamma*M/(1+(gamma-1)*M^2))/(c*k_p/alpha_p+n^2*c/r^2/k_0/alpha_p-gamma*M/(1+(gamma-1)*M^2));
            R_E_1=-(1-M^2)*(gamma-1)/(1+M^2*(gamma-1))*(1+R1);
            R_V_1=-c/k_0/r*(n*R1/r/alpha_p+n/r/alpha_n);
        else
            h = msgbox('This inlet boundary cannot be considered!');
        end
    end
end

% Outlet 
M=CI.TP.M_mean(1,end);
c=CI.TP.c_mean(1,end);
u=CI.TP.u_mean(1,end);
gamma=CI.TP.gamma(1,end);
% Defination of k and alpha in the duct
k_p=(M*omega-(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
k_n=(M*omega+(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
alpha_p=omega+u*k_p;
alpha_n=omega+u*k_n;
if strcmp(CI.setup.outlet,'open')==1
    R2=-1;
    R_E_2=0;
    R_V_2=0;
else if strcmp(CI.setup.outlet,'closed')==1 
        R2=-k_p/alpha_p/(k_n/alpha_n);
        R_E_2=0;
        R_V_2=0;
    else if strcmp(CI.setup.inlet,'choked')==1  % Stow and Dowling model 2002 -- for compact nozzle, l --> 0
            R2=-(c*k_p/alpha_p+0.5*(gamma-1)*M)/(c*k_n/alpha_n+0.5*(gamma-1)*M);
            R_E_2=-0.5*M/(c*k_n/alpha_n+0.5*(gamma-1)*M);
            R_V_2=n/(c*k_n/alpha_n+0.5*(gamma-1)*M);
            else
            h = msgbox('This outlet boundary cannot be considered!');
        end
    end
end
Te = Fcn_TF_entropy_convection(s); % Entropy dissipation/dispersion model
Tv = CI.setup.Tv;                  % Vorticity dissipation model
%% Index control
indexHA=0;

%% Initial wave amplitudes
A_0_n=1/R1;
A_0_p=1;
A_E_0=A_0_n*R_E_1;
A_V_0=A_0_n*R_V_1;
Array_LeftBD=[A_0_p, A_0_n, A_E_0, A_V_0].';
G = eye(4);
%% Transfer matrixes in each duct
% Calculate tranfer matrix from inlet to outlet
for ss = 1:CI.TP.numSection-1
    %% Calculate transfer matrixes for each duct and interface
    M=CI.TP.M_mean(1,ss);
    c=CI.TP.c_mean(1,ss);
    u=CI.TP.u_mean(1,ss);
    rho=CI.TP.rho_mean(1,ss);
    L_d=CI.setup.x(ss+1)-CI.setup.x(ss); % Duct length
    % Defination of k and alpha in the duct
    k_p=(M*omega-(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
    k_n=(M*omega+(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
    k_0=-omega/u;
    k_n0_p=-omega/(c+u);
    k_n0_m=-omega/(c-u);
    
    alpha_p=omega+u*k_p;
    alpha_n=omega+u*k_n;
    % Wave propagation matrix from left to right within the duct
    if CI.setup.InterfaceIndex(ss)==10 || CI.setup.InterfaceIndex(ss)==11  % Entropy and Vorticity dissipation/dispersion model
       P_l2r=[exp(1i*k_p*L_d), 0, 0, 0;
              0, exp(1i*k_n*L_d), 0, 0;
              0, 0, Te*exp(1i*k_0*L_d), 0;
              0, 0, 0,Tv*exp(1i*k_0*L_d)]; 
    else
        P_l2r=[exp(1i*k_p*L_d), 0, 0, 0;
              0, exp(1i*k_n*L_d), 0, 0;
              0, 0, exp(1i*k_0*L_d), 0;
              0, 0, 0,exp(1i*k_0*L_d)];
    end
    % Transfer matrix between wave amplitudes and flow perturbations in the duct
    T_w2f    = [1, 1, 0, 0;
                1/c^2, 1/c^2, -1/c^2, 0;
                -k_p/(rho*alpha_p), -k_n/(rho*alpha_n), 0, n/(rho*c);
                -n/(r*rho*alpha_p), -n/(r*rho*alpha_n), 0, -k_0*r/(rho*c)];
            
    % Inverse of T_w2f at the inlet of the next duct
      % After the interface
      M=CI.TP.M_mean(2,ss+1);
      c=CI.TP.c_mean(2,ss+1);
      u=CI.TP.u_mean(2,ss+1);
      rho=CI.TP.rho_mean(2,ss+1);
      k_p=(M*omega-(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
      k_n=(M*omega+(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
      k_0=-omega/u;
      alpha_p=omega+u*k_p;
      alpha_n=omega+u*k_n;
      T_w2f_inv_ai=[-(alpha_p*n^2 + alpha_p*k_0*k_n*r^2)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0,    -(alpha_n*alpha_p*k_0*r^2*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),                    -(alpha_n*alpha_p*n*r*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2);
                    (alpha_n*n^2 + alpha_n*k_0*k_p*r^2)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0,     (alpha_n*alpha_p*k_0*r^2*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),                     (alpha_n*alpha_p*n*r*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2);
                    1, -c^2,                                                                                                                                0,                                                                                                                                          0;
                    (c*k_n*n - c*k_p*n)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0, (rho*(alpha_n*c*n - alpha_p*c*n))/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2), (r*rho*(alpha_p*c*k_n - alpha_n*c*k_p))/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2)];  
      % After the flame
      switch CI.setup.InterfaceIndex(ss+1)
          case 0
              % Nothing happens
          case {10,11}
              M=CI.TP.M_mean(1,ss+1);
              c=CI.TP.c_mean(1,ss+1);
              u=CI.TP.u_mean(1,ss+1);
              rho=CI.TP.rho_mean(1,ss+1);
              % Defination of k and alpha in the duct
              k_p=(M*omega-(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
              k_n=(M*omega+(omega^2-(n^2*c^2/r^2)*(1-M^2))^(1/2))/c/(1-M^2);
              k_0=-omega/u;
              alpha_p=omega+u*k_p;
              alpha_n=omega+u*k_n;
              T_w2f_inv_af=[-(alpha_p*n^2 + alpha_p*k_0*k_n*r^2)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0,    -(alpha_n*alpha_p*k_0*r^2*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),                    -(alpha_n*alpha_p*n*r*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2);
                            (alpha_n*n^2 + alpha_n*k_0*k_p*r^2)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0,     (alpha_n*alpha_p*k_0*r^2*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),                     (alpha_n*alpha_p*n*r*rho)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2);
                            1, -c^2,                                                                                                                                0,                                                                                                                                          0;
                            (c*k_n*n - c*k_p*n)/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2),     0, (rho*(alpha_n*c*n - alpha_p*c*n))/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2), (r*rho*(alpha_p*c*k_n - alpha_n*c*k_p))/(alpha_n*n^2 - alpha_p*n^2 - alpha_p*k_0*k_n*r^2 + alpha_n*k_0*k_p*r^2)];  
      end
    %% Multiply transfer matrixes in downstream propagating direction
    DuctIndex=CI.setup.DuctIndex;
    DuctIndex=[DuctIndex,DuctIndex(end)];
    if DuctIndex(ss)==0 && DuctIndex(ss+1)==0  % Normal annular duct connections
        switch CI.setup.InterfaceIndex(ss+1)
            case 0
                Z       = T_w2f_inv_ai*CI.TPM.B2_inv{2,ss}*CI.TPM.B1{2,ss}*T_w2f*P_l2r;
            case 11
                indexHA = indexHA + 1;
                FTF     = Fcn_flame_model(s,indexHA);
                B2b     = zeros(4);
                B2b(4,3)= CI.TP.rho_mean(2,ss+1)*CI.TP.Q(indexHA)/CI.TP.mass(indexHA)*FTF;
                CI_TPM_B1_bf=CI.TPM.B1{1,ss}+B2b;
                Z       = T_w2f_inv_af*CI.TPM.B2_inv{1,ss}*CI_TPM_B1_bf*CI.TPM.B2_inv{2,ss}*CI.TPM.B1{2,ss}*T_w2f*P_l2r;
        end
    else if DuctIndex(ss)==0 && DuctIndex(ss+1)==1                          % From annular duct to premixers
            M_a2p=[1,0,0,0; 0,1,0,0; 0,0,0,1];
            Z    = CI.TPM.F_r_inv{2,ss}*CI.TPM.B2_inv{2,ss}*gamma_r_n*M_a2p*CI.TPM.B1{2,ss}*T_w2f*P_l2r;
        else if DuctIndex(ss)==1 && DuctIndex(ss+1)==0                      % From premixers back into annular duct
                P_l2r=[exp(1i*k_n0_p*L_d), 0, 0;
                       0, exp(-1i*k_n0_m*L_d), 0;
                       0, 0, exp(1i*k_0*L_d)];
                G_l  =CI.TPM.B1{2,ss}./gamma_r_n+CI.TPM.B1_1{2,ss};   
                M_p2a=[1,0,0; 0,1,0; 0,0,0; 0,0,1]; 
                %M_vortex_diss_test=[1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
                Z = T_w2f_inv_ai*CI.TPM.B2_inv{2,ss}*M_p2a*G_l*CI.TPM.F_l{2,ss}*P_l2r;
                else if DuctIndex(ss)==1 && DuctIndex(ss+1)==1                  % From premixers to premixers
                        P_l2r=[exp(1i*k_n0_p*L_d), 0, 0;
                            0, exp(-1i*k_n0_m*L_d), 0;
                            0, 0, exp(1i*k_0*L_d)];
                        Z = CI.TPM.F_r_inv{2,ss}*CI.TPM.B2_inv{2,ss}*CI.TPM.B1{2,ss}*CI.TPM.F_l{2,ss}*P_l2r;
                        % Other ductindex choices can be considered by adding other
                        % cases here
                    end
            end
        end
    end
    G=Z*G;
end
%% Outlet wave amplitudes
Array_RightBD   = G*Array_LeftBD;
A_p_N           = Array_RightBD(1);
A_n_N           = Array_RightBD(2);
A_E_N           = Array_RightBD(3);
A_V_N           = Array_RightBD(4);
F = ((R2*A_p_N + R_E_2*A_E_N+R_V_2*A_V_N) - A_n_N)/(abs(A_n_N)+abs(A_p_N)+abs(A_0_n)+abs(A_0_p));
%
%----------------------Entropy convection transfer function ---------------
%
function Te = Fcn_TF_entropy_convection(s)
global CI
tau     = CI.setup.tau_Cs;
k       = CI.setup.ke;
Te      = k*exp((tau*s)^2/4);

%
% ----------------------linear Flame transfer function --------------------
%         
function F = Fcn_flame_model(s,indexHA)
global CI
F = CI.setup.FM.a_f*exp(-CI.setup.FM.tau_f*s);
%
% -----------------------------end-----------------------------------------
