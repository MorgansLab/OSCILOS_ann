%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    OSCILOS-ann    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION OF COMBUSTION INSTABILITIES IN AN ANNULARCOMBUSTOR.
% THIS CODE INCLUDES BOTH A LINEAR AND A NONLINAER VERSIONS.
% Last update by Dong Yang, 02/12/2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add search paths
addpath(genpath('./'));  % Add the current path and all the subpaths to the MATLAB search path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%% Declare global variables
global CI

%% System setup                                        
CI=System_setup;                            % Edit this function in ./System_setup/ to change the system configurations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate mean flow parameters
CI=Fcn_calculation_mean_main(CI);           % Calculate mean flow profiles

% Calculate the transfer matrix from inlet to outlet
Fcn_PreProcessing;                          % Calculate transfer matrixes relating to only mean flow parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Frequency and growth rate ranges in which the thermoacoustic modes are calculated.
CI.EIG.Scan.FreqMin  =100;                  % The minimum frequency (Hz)
CI.EIG.Scan.FreqMax  =1000;                 % The maxmum frequency (Hz)
CI.EIG.Scan.GRMin    =-500;                 % The minimum growth rate (1/s)
CI.EIG.Scan.GRMax    =500;                  % The maximum growht rate (1/s)

%% Thermoacoustic mode calculation
% Choose a simulation model 
CI.CalStyle=2;                              % 1 denotes the linearly uncoupled model
                                            % 2 denotes the nonlinearly coupled model
switch CI.CalStyle
    case 1                                  
    % The linearly uncoupled model. Set the scan domain for the frequency and growth rate
    CI.setup.n=1;                           % Fixed circumferential wave number (assume that there is no modal coupling)
    CI.EIG.Scan.FreqNum  =10;               % The number for initial frequency guess within the given frequeny range
    CI.EIG.Scan.GRNum    =10;               % The number of initial growth rate guess within the given growth rate range
    case 2
    % The nonlinearly coupled model
        CI.setup.N=3;                                               % Truncation number of circumferential mode number
        CI.fixed_growthrate  =0;                                    % This is the fixed growth rate of the targeted mode (1/s)    
        CI.f_iniguess        =481;                                  % Initial guess of the angular frequency (Hz)
        CI.lambda_iniguess   =[0,0,3.2,0,0,0,6.9,3.2,0,0,0,0,0];  % Initial guess of the modal amplitudes (Pa)
    otherwise
        h = msgbox('This CalStyle is not considered currently!');
end
assignin('base','CI',CI)
Eigenmode         = Fcn_calculation_eigenmode;  
CI.Eigenmode.modes= Eigenmode;
CI.Eigenmode.GR   = real(Eigenmode);        % Growth rates of the modes
CI.Eigenmode.Freq = imag(Eigenmode)/2/pi;   % Frequencies of the modes    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Result
Figure_num            =1;                   % Define a figure number
switch CI.CalStyle
    case 1                                  % For the linearly uncoupled model
        % Contour plot
        Fcn_contour_plot(Figure_num);       % Contour plot over the scan frequency and growth rate range
        % Modeshape
        Figure_num    =Figure_num+1;
        Mode_num      =1;                   % The mode number to plot
        Angle2plot    =90;                  % Circumferential angle to plot (degree)
        Fcn_osc_modeshape_linear_n_V0(Figure_num,Mode_num,Angle2plot);
    case 2                                  % For the nonlinearly coupled case
        % Modeshape plot (along the longitudinal dirention, at a given fix angel)
        Figure_num    =Figure_num+1;        % Define a figure number
        Angle2plot    =0;                   % Circumferential angle to plot (degree)
        n_2plot       =[-1:1:1];            % Circumferential modal components to plot
        Fcn_osc_modeshape_nonlinear_V0(Figure_num,Angle2plot,n_2plot);
        % Perturbations just at the burners' outlets (before the flames), and the flame
        % heat release oscillations
        Figure_num    =Figure_num+1;
        Fcn_osc_perturbations_at_burner_outlets(Figure_num); 
        % Track the evolution of the nonlinearly coupled eigenmode (such as from an unstable mode solution to a limit cycle)
        Stop_growth_rate=100;               % The evolution starts from the mode with a growth rate =CI.fixed_growthrate, and try to evolute towards a mode with this Stop_growth_rate
        Iteration_steps=20;
        S_nonlinear = Fcn_calculation_tracking_one_nonlinear_mode(Iteration_steps, Stop_growth_rate);
        Figure_num=Figure_num+1;
        Fcn_nonlinear_mode_evolution(S_nonlinear, Figure_num) % Plot the mode evolution
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove search paths
rmpath(genpath('./'));  % Remove the current path and all the subpaths to the MATLAB search path
% -----------------------------end-----------------------------------------