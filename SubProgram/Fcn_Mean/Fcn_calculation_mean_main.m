function CI=Fcn_calculation_mean_main(CI)
% This function is used to calculate the mean Thermal Properties in every section
% first created: 05/10/2016
% last modified: 02/12/2018
% ----------------------------
CI.TP.numSection=length(CI.setup.InterfaceIndex);
% -----------calculate the properties in different sections----------------
for ss = 1:CI.TP.numSection-1 
    CI.TP.Theta(ss)  = CI.setup.Area(ss+1)./CI.setup.Area(ss);       % Surface area ratio S2/S1
end
% --------------------------------
% Inlet mean flow parameters
CI.TP.p_mean(1:2,1)  =CI.setup.p1;
CI.TP.rho_mean(1:2,1)=CI.setup.rho1;
CI.TP.u_mean(1:2,1)  =CI.setup.u1;
CI.TP.T_mean(1:2,1)  =CI.setup.T1;
CI.TP.gamma(1:2,1)   =CI.setup.gamma_bf;
CI.TP.Cp(1:2,1)      =CI.setup.Cp;
CI.TP.c_mean(1:2,1)  =CI.setup.c1;
CI.TP.M_mean(1:2,1)  =CI.setup.M1;

indexHA_num=0;         % set the initial value to zero and it will be increased by 1 after every HA interface
for ss = 1:CI.TP.numSection-1 
    % In every interface, the changes are splitted to two steps:
    % 1. cross sectional surface area change
    % 2. Heat addition or .....
    % --------------step 1-------------------------
    %
    [   CI.TP.p_mean(1:2,ss+1),...
        CI.TP.rho_mean(1:2,ss+1),...
        CI.TP.u_mean(1:2,ss+1) ]... 
      = Fcn_calculation_mean_area_change(CI.TP.p_mean(1,ss),...
                                         CI.TP.rho_mean(1,ss),...
                                         CI.TP.u_mean(1,ss),...
                                         CI.TP.Theta(ss),...
                                         CI.TP.gamma(1,ss));
    % ----------
    CI.TP.gamma(1:2,ss+1)   = CI.TP.gamma(1,ss);
    CI.TP.Cp(1:2,ss+1)      = CI.TP.Cp(1,ss);
    CI.TP.T_mean(1:2,ss+1)  = CI.TP.gamma(1,ss+1)/(CI.TP.gamma(1,ss+1)-1)...
                             *CI.TP.p_mean(1,ss+1)./(CI.TP.Cp(1,ss+1).*CI.TP.rho_mean(1,ss+1));
    CI.TP.c_mean(1:2,ss+1)  = ((CI.TP.gamma(1,ss+1) - 1).*CI.TP.Cp(1,ss+1).*CI.TP.T_mean(1,ss+1)).^0.5;
    CI.TP.M_mean(1:2,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
    %
    % --------------step 2-------------------------
    %
    switch CI.setup.InterfaceIndex(ss+1)
        case 0
            % in case 0, no changes
        case {11}                                       % with HA
            indexHA_num = indexHA_num + 1;              % this number is increased by 1
            %Calculate temperature, gamma and R_air after the flame
            CI.TP.T_mean(1,ss+1)  = CI.setup.T_af;
            CI.TP.gamma(1,ss+1)   = CI.setup.gamma_af;
            CI.TP.Cp(1,ss+1)      = CI.setup.gamma_af/(CI.setup.gamma_af-1)*CI.setup.R_air;
                    % ---------then, use the resolved temperature, gamma and the mean
                    % properties after the area changes to calculate the mean
                    % properties after HA ----------------------------------------
                    [   CI.TP.p_mean(1,ss+1),...               
                        CI.TP.rho_mean(1,ss+1),...
                        CI.TP.u_mean(1,ss+1)] = ...
                        Fcn_calculation_mean_across_flame(CI.TP.p_mean(2,ss+1),...
                                                          CI.TP.rho_mean(2,ss+1),...
                                                          CI.TP.u_mean(2,ss+1),...
                                                          CI.setup.R_air,...
                                                          CI.TP.T_mean(1,ss+1));
                    CI.TP.c_mean(1,ss+1)  = sqrt(CI.TP.gamma(1,ss+1)*CI.setup.R_air*CI.TP.T_mean(1,ss+1));
                    CI.TP.M_mean(1,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
                    
                    CI.TP.DeltaHr(indexHA_num) = CI.TP.Cp(1,ss+1).*(CI.TP.T_mean(1,ss+1) - CI.TP.T_mean(2,ss+1))...
                        + 0.5*(CI.TP.u_mean(1,ss+1).^2 - CI.TP.u_mean(1,ss).^2);              
                    CI.TP.mass(indexHA_num)    = CI.TP.rho_mean(2,ss+1).*CI.TP.u_mean(2,ss+1).*CI.setup.Area(ss+1); % mass flow rate before HA
                    CI.TP.Q(indexHA_num)       = CI.TP.DeltaHr(indexHA_num).*CI.TP.mass(indexHA_num);               % heat release rate             
    end
end
% -------------------------------end --------------------------------------