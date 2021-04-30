function F = Fcn_Q_n0(s,lambda,ss_premixer)

global CI

N_all=2*CI.setup.N+1;
n_all=-CI.setup.N:1:CI.setup.N;

D=CI.setup.Premixer_Num;
phi_d_all_coe=0:1/D:(D-1)/D;
phi_d_all=2*pi*phi_d_all_coe;


[p_3_n,u_3_n]=Fcn_u_3_n(s,ss_premixer); % To calculate velocity oscillations before the flames.
u_3_n=u_3_n.*lambda;            % This is where lambda comes in


u_mean_d=CI.TP.u_mean(1,ss_premixer);   % This mean velocity is assumed to be the same for all premixers. This is inside the premixers 
%u_mean_d=CI.TP.u_mean(1,ss_premixer+1);   % This mean velocity is assumed to be the same for all premixers. This is at the annular duct just after the premixers 

indexHA=1;                              % Only one flame interface is considered in the present version of OSCILOS_Annular
Q_mean_d=CI.TP.Q(indexHA)/D;            % This mean heat release rate is assumed to be the same for all premixers.


for n_num=1:1:N_all
    n0=n_all(n_num);

for d=1:1:D
    A_u_d(d) = sum(u_3_n.*exp(1i.*n_all.*phi_d_all(d)))/u_mean_d;
    Flame_model_d(d)=Fcn_flame_model(s,A_u_d(d),indexHA);
    % Q_d
    Q_d(d)   = Q_mean_d*Flame_model_d(d)*A_u_d(d);
    % Q_n0_d
    Q_n0_d(d)= Q_mean_d/u_mean_d*Flame_model_d(d)*sum(u_3_n.*exp(1i.*(n_all-n0).*phi_d_all(d)));
end
Q_n0(n_num)=sum(Q_n0_d);

end

F=Q_n0;  
%       