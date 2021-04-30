function [p_mean2,rho_mean2,u_mean2] =...
    Fcn_calculation_mean_area_change(p_mean1,rho_mean1,u_mean1,Theta,gamma)
% This function is used to calculate the mean thermoproperties changing
% across the interface between two sections with different section surface
% area. Heat addition is not accounted for, gamma is constant.
% Theta=S2/S1;
% 
MF.p1       = p_mean1;
MF.rho1     = rho_mean1;
MF.u1       = u_mean1;
MF.Theta    = Theta;
MF.gamma    = gamma;
[p_mean2, rho_mean2, u_mean2] = FcnSolve(@myfun_p_rho_u, MF);
%
function [p_mean2, rho_mean2, u_mean2] = FcnSolve(myfun, MF)
%
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
x0      = [1,1,1];  
F       = @(x)myfun(x, MF);
iternum=0;
exitflag=0;
while exitflag>4 || exitflag<1 && iternum<20
    [x2,fval,exitflag]     = fsolve(F, x0, options);   % solve equation
    iternum=iternum+1;
    x0=x2;
end
if iternum>19
    h = msgbox('The mean flow calculation is not converged!'); 
end
p_mean2     = x2(1)*MF.p1;               
rho_mean2   = x2(2)*MF.rho1;
u_mean2     = x2(3)*MF.u1;
%
function F = myfun_p_rho_u(x, MF)
%
p1      = MF.p1;
rho1    = MF.rho1;
u1      = MF.u1;
Theta   = MF.Theta;
gamma   = MF.gamma;
%
p2      = x(1)*p1;
rho2    = x(2)*rho1;
u2      = x(3)*u1;
F(1) = Theta*rho2*u2 - rho1*u1;
%
F(3) = gamma/(gamma-1)*(p2/rho2 - p1/rho1) + 0.5*(u2^2 - u1^2);
%
if Theta>=1              % Area increasing
    F(2) = p2+rho2*u2^2 - (p1 + rho1*u1^2/Theta);
elseif Theta<1          % Area decreasing
    F(2) = p2*rho1^gamma - p1*rho2^gamma;
end
%
%----------------------------------end-------------------------------------



%
% [p_mean2, rho_mean2, u_mean2] = FcnSolve(MF);
% %
% function [p_mean2, rho_mean2, u_mean2] = FcnSolve(MF)
% %
% if MF.Theta>=1
%     % Area increase: Coeff_2_34*u_mean2^2+Coeff_1_34*u_mean2+Coeff_0_34=0
%     Coeff_2_34=(1+MF.gamma)/2/(1-MF.gamma);
%     Coeff_1_34=MF.gamma/(MF.gamma-1)*(MF.p1*MF.Theta/(MF.rho1*MF.u1)+MF.u1);
%     Coeff_0_34=-(MF.gamma/(MF.gamma-1)*MF.p1/MF.rho1)-0.5*MF.u1^2;
%     syms t4 positive;
%     u_mean2=solve(Coeff_2_34*t4^2+Coeff_1_34*t4+Coeff_0_34,'PrincipalValue', true, 'Real',true); % Average flow velocity in the short chamber before premixers
%     u_mean2=double(u_mean2);
%     rho_mean2=MF.rho1*MF.u1/MF.Theta/u_mean2;                                             % Mean density in the short chamber before premixers
%     p_mean2=MF.p1+MF.rho1*MF.u1^2/MF.Theta-rho_mean2*u_mean2^2;                           % Mean pressure in the short chamber before premixers
% else % Area decrease
%     Coeff_2_12=0.5;
%     Coeff_1_12=MF.gamma/(MF.gamma-1)*MF.p1/MF.rho1*(MF.Theta/MF.u1)^(1-MF.gamma);
%     Coeff_0_12=-MF.gamma/(MF.gamma-1)*MF.p1/MF.rho1-0.5*MF.u1^2;
%     ufun12=@(t2,Coeff_2_12, Coeff_1_12, Coeff_0_12, gamma) Coeff_2_12*t2^2+Coeff_1_12*t2^(1-gamma)+Coeff_0_12;
%     t2=[0.1,500]; % This range may need to be valified for a specific mean flow profile
%     myfun12=@(t2) ufun12(t2,Coeff_2_12, Coeff_1_12, Coeff_0_12, MF.gamma);
%     u_mean2=fzero(myfun12,t2);
%     rho_mean2=MF.rho1*MF.u1/MF.Theta/u_mean2;                                             % Mean density in the inlet duct
%     p_mean2=(rho_mean2/MF.rho1)^(MF.gamma)*MF.p1;                                         % Mean pressure in the inlet duct
% end
% %
% %----------------------------------end-------------------------------------
