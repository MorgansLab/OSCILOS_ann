% ----------------------------Calculate eigenvalue-------------------------
function [E,initial_guess_new] = Fcn_calculation_eigenmode
% This function is used to calculate the eigenvalue 
global CI
%---------------------------------
GRMin       = CI.EIG.Scan.GRMin;   
GRMax       = CI.EIG.Scan.GRMax;
FreqMin     = CI.EIG.Scan.FreqMin;   
FreqMax     = CI.EIG.Scan.FreqMax;

switch CI.CalStyle
    case 1
        GRNum       = CI.EIG.Scan.GRNum;
        GRSp        = linspace(GRMin,   GRMax,      GRNum);         % sample number of growth rate
    case 2
        GRSp        = linspace(0, 0, 1);
end
switch CI.CalStyle
    case 1
        FreqNum     = CI.EIG.Scan.FreqNum;
        FreqSp      = linspace(FreqMin, FreqMax,    FreqNum);       % sample number of scan frequency
    case 2
        FreqSp      = linspace(0, 0, 1);
end
%---------------------------------
eigen_num=1;
h = waitbar(0,'Please wait...');
for kk = 1:length(FreqSp) 
    omega = 2*pi*FreqSp(kk);
    for ss = 1:length(GRSp)
        GR = GRSp(ss);
        if GR == 0
            GR = 1;
        end
        if omega == 0
            omega =10;
        end
        s0 = GR+1i*omega;                                                          % initial value
        switch CI.CalStyle
            case 1                                                                 % Continuous flame model -- can be linear or nonlinear, but modal interactions are not considered
                options = optimset('Display','off','TolFun',1e-25,'TolX',1e-10);   % the calculation results by fsolve will not be shown in the workspace
                [x,fval,exitflag] = fsolve(@Fcn_DetEqn_Linear,s0,options); 
                x(1)=x(1);
            case 2                                                                 % Discreted flame model -- modal interactions are considered; D>2N+1 where D is the flame number and N is the acoustic modal truncation number          
                omega =CI.f_iniguess*2*pi;
                lambda=CI.lambda_iniguess;
                Ini_guess(1)=omega;
                Ini_guess(2:(1+CI.setup.N))=lambda(1:CI.setup.N);
                Ini_guess((2+CI.setup.N):(1+2*CI.setup.N))=lambda((1+CI.setup.N):2*CI.setup.N);
                Ini_guess(2+2*CI.setup.N)=lambda(1+2*CI.setup.N);
                Ini_guess((3+2*CI.setup.N):(2+3*CI.setup.N))=lambda((2*CI.setup.N+2):(3*CI.setup.N+1));
                Ini_guess((3+3*CI.setup.N):(2+4*CI.setup.N))=lambda((3*CI.setup.N+2):(4*CI.setup.N+1));

                options = optimset('MaxIter',100,'MaxFunEvals',10000,'Display','iter','TolFun',1e-15,'TolX',1e-7);        % the calculation results by fsolve will not be shown in the workspace
                [x,fval,exitflag,output,jacobian] = fsolve(@Fcn_DetEqn_NonLinear,Ini_guess,options);
                x=x;
                fval=fval;
                exitflag=exitflag;
                jacobian=jacobian;
        end
         if exitflag>0
            if CI.CalStyle==2  
                EIG.eigenvalue(eigen_num) = 1i*x(1)+CI.fixed_growthrate; % To ensure the same frequency form as that is case 1
                EIG.nonlinear_solution{eigen_num}=x;
                initial_guess_new=x(2:end);
            else                                                                  % This can only be CI.CalStyle==1 in the present version of code
                EIG.eigenvalue(eigen_num) = x(1);
            end
            EIG.eigenvalue_prec(eigen_num) = round(EIG.eigenvalue(eigen_num).*100)./100;      % this is used to set the precision
            eigen_num = eigen_num+1;
         end
    end
     waitbar(kk / length(FreqSp))  
end
close(h)
[b,m,n] = unique(EIG.eigenvalue_prec);                      % unique function is used to get rid of the same value
EIG.eigenvalue_unique = EIG.eigenvalue(m);                  % this is the eigenvalue
if CI.CalStyle==2
    EIG.nonlinear_solution_unique=EIG.nonlinear_solution(m);
else
end
EIG.eigenvalue_unique_prec = EIG.eigenvalue_prec(m);
%---------------------------------
% this is used to get rid of the zero frequency value
s_null = [];
cal = 0;
for ss = 1:length(EIG.eigenvalue_unique)
    if(abs(imag(EIG.eigenvalue_unique(ss)))<1)
        cal = cal+1;
        s_null(cal) = ss;
    end
end
EIG.eigenvalue_unique(s_null) = [];
if CI.CalStyle==2
    EIG.nonlinear_solution_unique(s_null)= [];
else 
end
%---------------------------------
% the previous processing is still not enough to get rid of the same
% value,such as 999.9 and 1000.1 corresponding to the function `floor'
% such as 995.1 and 994.9 corresponding to the function `round'
if length(EIG.eigenvalue_unique)>1
    cal = 0;
    EIG.eigenvalue_unique = sort(EIG.eigenvalue_unique);
    EIG.eigenvalue_unique_diff = diff(EIG.eigenvalue_unique);
    EIG.index_NULL = [];
    for kk = 1:length(EIG.eigenvalue_unique_diff)
        if(abs(EIG.eigenvalue_unique_diff(kk))<10)
            cal = cal+1;
            EIG.index_NULL(cal) = kk;
        end
    end
    EIG.eigenvalue_unique(EIG.index_NULL) = []; % this is the eigenvalue we want
    if CI.CalStyle==2
        EIG.nonlinear_solution_unique(EIG.index_NULL) = [];
    else
    end
else 
end
EIG.growthrate_limit    = [GRMin    GRMax];
EIG.Omega_limit         = [FreqMin  FreqMax].*2.*pi;
%-------------------------------------
s_null=[];
cal=0;
for ss=1:length(EIG.eigenvalue_unique)
    if(real(EIG.eigenvalue_unique(ss))<EIG.growthrate_limit(1)||real(EIG.eigenvalue_unique(ss))>EIG.growthrate_limit(2)||abs(imag(EIG.eigenvalue_unique(ss)))<EIG.Omega_limit(1)||abs(imag(EIG.eigenvalue_unique(ss)))>EIG.Omega_limit(2) )
        cal=cal+1;
        s_null(cal)=ss;
    end
end
EIG.eigenvalue_unique(s_null)=[];        % Do not use semicolon if you want to show the final value
if CI.CalStyle==2
    EIG.nonlinear_solution_unique(s_null) = [];
else
end
E = EIG.eigenvalue_unique; 
if CI.CalStyle==2
    CI.EIG.nonlinear_solution_unique=EIG.nonlinear_solution_unique;
    assignin('base','CI',CI)
else 
end
clear EIG