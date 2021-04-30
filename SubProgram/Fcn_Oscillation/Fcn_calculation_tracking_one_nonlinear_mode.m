% ----------------------------Calculate eigenvalue evolution-------------------------
function S_nonlinear = Fcn_calculation_tracking_one_nonlinear_mode(Iteration_steps, Stop_growth_rate)
% This function is used to track the given nonlinear eigenvalue 
global CI
%% Give one of the mode we want to tract -- this comes from previous nonlinear calculation (needs to see the change of growth rate, otherwise the iteration may not work)
s_nonlinear_0=CI.EIG.nonlinear_solution_unique{1};
growth_rate_0=CI.fixed_growthrate;
growth_rate_step=(Stop_growth_rate-growth_rate_0)/Iteration_steps;

s_nonlinear=s_nonlinear_0;
h = waitbar(0,'Please wait...');
for N_iter=1:1:(Iteration_steps+1)
    growth_rate=growth_rate_0+growth_rate_step*(N_iter-1);
    CI.fixed_growthrate=growth_rate;
    assignin('base','CI',CI);
    
    omega =s_nonlinear(1);
    s_lambda=s_nonlinear(2:(4*CI.setup.N+2));
    Ini_guess(1)=omega;
    Ini_guess(2:(1+CI.setup.N))=s_lambda(1:CI.setup.N);
    Ini_guess((2+CI.setup.N):(1+2*CI.setup.N))=s_lambda((1+CI.setup.N):2*CI.setup.N);
    Ini_guess(2+2*CI.setup.N)=s_lambda(1+2*CI.setup.N);
    Ini_guess((3+2*CI.setup.N):(2+3*CI.setup.N))=s_lambda((2*CI.setup.N+2):(3*CI.setup.N+1));
    Ini_guess((3+3*CI.setup.N):(2+4*CI.setup.N))=s_lambda((3*CI.setup.N+2):(4*CI.setup.N+1));
    
    options = optimset('MaxIter',400,'MaxFunEvals',10000,'Display','off','TolFun',1e-5,'TolX',1e-4);        % the calculation results by fsolve will not be shown in the workspace
    [x,fval,exitflag,output,jacobian] = fsolve(@Fcn_DetEqn_NonLinear,Ini_guess,options);
    exitflag=exitflag;
    if exitflag>0
        s_nonlinear=x;
        S_nonlinear(N_iter,1)=growth_rate;
        S_nonlinear(N_iter,2:(1+length(s_nonlinear)))=s_nonlinear;
    else
        h = msgbox('The growth rate step is too big!');
        break
    end
    waitbar(N_iter / (Iteration_steps+1))  
end  
close(h)