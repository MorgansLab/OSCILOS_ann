      
function F = Fcn_flame_model(s,A_u_d,indexHA)
% ----------------------Flame transfer function --------------------
%   
global CI

switch CI.setup.FM.NStyle
    case 1
        %% Stow & Dowling's (ASME2009) nonlinear flame model
        if abs(A_u_d)>CI.setup.FM.SD.phaselim
            tau_f=CI.setup.FM.tau_f*(1+CI.setup.FM.SD.beta*(abs(A_u_d)-CI.setup.FM.SD.phaselim));          % Nonconstant time delay
        else
            tau_f=CI.setup.FM.tau_f;
        end
        F        = CI.setup.FM.a_f*exp(-tau_f*s);
        
        beta     =abs(CI.setup.FM.a_f*A_u_d)/CI.setup.FM.SD.alpha;
        if beta<=1
            F=F;
        else
            F=F*(1-2*acos(1/beta)/pi+2*sqrt(1-1/beta^2)/(pi*beta));
        end
    case 2
        %% Li & Morgans' (JSV2015) nonlinear flame model
        fun    =@(x) 1./(1+(x+CI.setup.FM.LM.alpha).^CI.setup.FM.LM.beta);
        L      =integral(fun,0,abs(A_u_d))/abs(A_u_d);
        tau_f  =CI.setup.FM.tau_f+CI.setup.FM.LM.tau_fN*(1-L);
        omega_c=2*pi*CI.setup.FM.LM.fc;
        F      = L*(CI.setup.FM.a_f*omega_c)/(s+omega_c)*exp(-tau_f*s);
end