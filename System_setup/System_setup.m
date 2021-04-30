% --- Creat CI.S_setup. to give the system setup
function CI=System_setup(varargin)
% This program is used to setup all the input parameters for the system
%
% ---------------------------------
%% Combustor geometry
R_m   = 175/1000;                        % Mean radius of the annular duct (mm) -- assumed the same value for all ducts. 

L_ori = [70 14 6 0 241]./1000;           % Original longitudinal lengthes (mm)      
L     = [70 14+5 6+3 0 241]./1000;       % Lengthes with small endlength corrections (mm)
Area  = [0.076969 0.013685 0.0044736 0.054978 0.054978]; % Cross-sectional areas of ducts,(m^2) -- Areas for the premix ducts sections are the total area of all premix ducts
CI.setup.DuctIndex     = [0 1 1 0 0];    % Set and index for each duct.
                                         % "0" denotes normal annular straight duct
                                         % "1" denotes premix ducts (cylinders)
CI.setup.Premixer_Num  = 16;             % The number of premix ducts/burners
CI.setup.InterfaceIndex= [0 0 0 0 11 0]; % Set an index for each interface, including the inlet and outlet
                                         % "0" denotes normal area expansion or contraction; 
                                         % "11" denotes heat additions (both mean and perturbation) 
CI.setup.x_ori= [0,cumsum(L_ori)];                                               
CI.setup.x    = [0,cumsum(L)];           % Genearate axial location data
CI.setup.Area = [Area,Area(end)];        % Generate  cross-sectional area data. Each value denotes the cross-sectional area just after the corresponding interface
CI.setup.R_m  = R_m;                     % Generate mean radius data

%% Boundary condition
CI.setup.inlet='closed';                 % Inlet boundary condition: chose from 'open', 'closed', and 'choked'
CI.setup.outlet='open';                  % Outlet boundary condition: chose from 'open', 'closed', and 'choked'

%% Thermodynamic data
CI.setup.R_air=287.1;                    % Ideal gas constant for air. Use the same value throughout the system (J/kg/K)

%% Inlet mean flow
CI.setup.T1=300;                         % Inlet mean temperature (K)
CI.setup.p1=1.0133E5;                    % Inlet mean pressure (Pa)
CI.setup.massflow1=8.478E-3;             % Inlet mean mass flow rate (kg/s)
% Calculate the rest inlet mean flow parameters
[~,~,~,Cp_bf]= Fcn_calculation_c_q_air(CI.setup.T1);
CI.setup.gamma_bf   = Cp_bf./(Cp_bf - CI.setup.R_air);              % Heat capacity ratio before flame (J/kg/K)
CI.setup.rho1=CI.setup.p1/(CI.setup.R_air*CI.setup.T1);             % Inlet mean density (kg/m^3)
CI.setup.c1=sqrt(CI.setup.gamma_bf*CI.setup.R_air*CI.setup.T1);     % Inlet mean sound speed (m/s)
CI.setup.u1=CI.setup.massflow1/CI.setup.rho1/CI.setup.Area(1);      % Inlet mean velocity (m/s)
CI.setup.M1=CI.setup.u1/CI.setup.c1;                                % Inlet mean Mach number
CI.setup.Cp=CI.setup.gamma_bf/(CI.setup.gamma_bf-1)*CI.setup.R_air; % Inlet heat capacitiy at constant pressure (J/kg/K)
    
%% Mean heat addition
CI.setup.T_af=1500;                                                 % Mean temperature after the flame (K)
[~,~,~,Cp_af]= Fcn_calculation_c_q_air(CI.setup.T_af);              % Heat capacitiy at constant pressure (J/kg/K)
CI.setup.gamma_af   = Cp_af./(Cp_af - CI.setup.R_air);              % Heat capacity ratio after flame
%% Oscillating heat addition (flame model)                                                 
% Linear n-tau model
CI.setup.FM.a_f=1;                       % Flame interaction index
CI.setup.FM.tau_f=1.0*10^(-3);           % Time delay between heat release rate and incident velocity, s 
                                         % Typical n-tau model is used here. Q'/Q_average=u'/u_average*(a_f*exp(-tau_f*s))                                             
                                               
% Nonlinear flame model -- this is needed when a nonlinear calculation is
% performed.
CI.setup.FM.NStyle       =2;             % "1" means using Stow and Dowling's (ASME2009) nonlinear model
                                         % "2" means using Li and Morgans (JSV2015) nonlinear model 
% Stow and Dowling's (ASME2009) nonlinear flame model incorporated with an amplitude
% dependent time delay
CI.setup.FM.SD.alpha    =0.5;            % Sturation constant in Stow's (ASME2009) nonlinear flame model
CI.setup.FM.SD.beta     =0.00000;        % Linear variation coefficient of the time lag; 0 means no amplitude dependance
CI.setup.FM.SD.phaselim =0.5;
% Li and Morgans' (JSV2015) nonlinear flame model with a first-order low pass filter
CI.setup.FM.LM.alpha  =0.9;              % This is the alpha in Li and Morgans' model
CI.setup.FM.LM.beta   =40;               % This is the beta in Li and Morgans' model
CI.setup.FM.LM.tau_fN =0.1*10^(-3);      % This is the coefficient of the time delay dependence on L, \tau_f^N, in Li and Morgans' model
CI.setup.FM.LM.fc     =100000;           % This is the cut-off frequency of the first-order low pass filter 
%% Entropy and vorticity propagation model after the flame
CI.setup.ke=0;                           % Entropy dissipation model; 1 means no dissipation and 0 totally dissipated
CI.setup.tau_Cs=0;                       % Entropy dispersion model; 0 means no dispersion
CI.setup.Tv=0;                           % Vorticity dissipation model; 1 means no dissipation and 0 totally dissipated
end                                     
% -----------------------------end-----------------------------------------