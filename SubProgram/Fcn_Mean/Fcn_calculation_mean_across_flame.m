function [p_mean2,rho_mean2,u_mean2]...
        = Fcn_calculation_mean_across_flame(p_mean1,rho_mean1,u_mean1,Rg2,T_mean2)
% This function is used to calculate the mean thermoproperties changing
% across the flame front which is considered as compact
% The calculation  consists two steps
% calculate the properties after the flame, herein, considere there is
% section areas are constant which equals to that of section 2
%-----------------------------step  ---------------------------------------
global MF
MF.x1   = [p_mean1,rho_mean1,u_mean1];
MF.Rg2  = Rg2;
MF.T2   = T_mean2;
%--------------------------------------------------------------------------
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
x0 = MF.x1;  % Make a starting guess at the solution
%
iternum=0;
exitflag=0;
while exitflag>4 || exitflag<1 && iternum<100
    [x2,fval,exitflag]     = fsolve(@myfun,x0,options); % Call solver
    iternum=iternum+1;
    x0=x2;
end
% if iternum>99
%     h = msgbox('The mean flow calculation is not converged!'); 
% end
p_mean2     = x2(1);               
rho_mean2   = x2(2);
u_mean2     = x2(3);
clear MF

function F = myfun(x2)
global MF
% p2    = x2(1);
% rho2  = x2(2);
% u2    = x2(3);
x1          = MF.x1;
Rg2         = MF.Rg2;
T2          = MF.T2;
F(1) = x2(2)*x2(3) - x1(2)*x1(3);                               % mass conservation equation
F(2) = (x2(1)+x2(2)*x2(3)^2) - (x1(1)+x1(2)*x1(3)^2);           % momentum conservation equation
F(3) = x2(1) - x2(2)*Rg2*T2;                                    % gas law equation
% ---------------------------------end-------------------------------------