function [W_amp_along_axial, P_amp_along_axial, x_discr] = Fcn_osc_modeshape_linear_n_V0(Figure_num,Mode_num,Angle2plot)
% This function only accounts for linear cases
global CI
s=CI.Eigenmode.modes(Mode_num);
omega=-1i*s;       % This depends on the defination of s.

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
else if strcmp(CI.setup.inlet,'closed')==1 % This may not be strictly accurate and needs to be validated latter.
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
else if strcmp(CI.setup.outlet,'closed')==1  % This is not strictly accurate and needs to be improved latter.
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
A_0_n=1;
A_0_p=A_0_n*R1;
A_E_0=A_0_n*R_E_1;
A_V_0=A_0_n*R_V_1;
Array_LeftBD=[A_0_p, A_0_n, A_E_0, A_V_0].';
Array_LeftBD=Array_LeftBD*exp(1i*n*Angle2plot/180*pi);
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
   
            
    %% Begin discretisation of each duct
    DuctIndex=CI.setup.DuctIndex;
    DuctIndex=[DuctIndex,DuctIndex(end)];
    x_all=CI.setup.x;
    x_diff=diff(x_all);                                                 % Section length
    Sec_dis_num=50;
    if x_diff(ss)==0
        x_section_discr=zeros(1,Sec_dis_num+1);
    else
        x_section_discr=0:x_diff(ss)/Sec_dis_num:x_diff(ss);
    end
    
    for sec_dis_num=1:1:length(x_section_discr)
        if DuctIndex(ss)==0                                             % This is normal annular duct
            P_l2x=[exp(1i*k_p*x_section_discr(sec_dis_num)), 0, 0, 0;
                0, exp(1i*k_n*x_section_discr(sec_dis_num)), 0, 0;
                0, 0, exp(1i*k_0*x_section_discr(sec_dis_num)), 0;
                0, 0, 0,exp(1i*k_0*x_section_discr(sec_dis_num))];
            T_w2f_no_V=T_w2f*[1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 0];
            W_amp_along_axial{sec_dis_num,ss}   = P_l2x*G*Array_LeftBD;
            P_amp_along_axial{sec_dis_num,ss}   = T_w2f_no_V*W_amp_along_axial{sec_dis_num,ss}; % Neglect the vortex effects
            
        else if DuctIndex(ss)==1
                P_l2x=[exp(1i*k_n0_p*x_section_discr(sec_dis_num)), 0, 0;
                    0, exp(-1i*k_n0_m*x_section_discr(sec_dis_num)), 0;
                    0, 0, exp(1i*k_0*x_section_discr(sec_dis_num))];
                T_w2f_no_V=CI.TPM.F_l{2,ss};
                W_amp_along_axial{sec_dis_num,ss}   = P_l2x*G*Array_LeftBD;
                P_amp_along_axial{sec_dis_num,ss}   = T_w2f_no_V*W_amp_along_axial{sec_dis_num,ss}; % There is no vortex in this section
            else
                % Left for other possible duct options
            end
        end
        x_discr{sec_dis_num,ss}             = CI.setup.x(ss)+x_section_discr(sec_dis_num);
    end
     %% End discretisation of each duct
                  
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
    else if DuctIndex(ss)==0 && DuctIndex(ss+1)==1                          %#ok<SEPEX> % From annular duct to premixers
            M_a2p=[1,0,0,0; 0,1,0,0; 0,0,0,1];
            %Z1  = CI.TPM.B1{2,ss}*T_w2f*P_l2r
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
                    else
                    end
            end
        end
    end      
    G=Z*G;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Modeshape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure(Figure_num);
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(3).*(1/100) scrsz(4).*(1/20) scrsz(3)*5/7 scrsz(4).*(6/7)])
set(h,'name','Mode analysis for a choked-choked straight duct','numbertitle','off');
%************
hAxes(1)=axes('Unit','pixels','position',[100 100 500 300]);
hold on
p_max=0;
for section_plot=1:1:CI.TP.numSection-1
    for sec_dis_num=1:1:length(x_section_discr)
        x_plot(sec_dis_num)=x_discr{sec_dis_num,section_plot}(1);
        p_plot(sec_dis_num)=abs(P_amp_along_axial{sec_dis_num,section_plot}(1));
        u_plot(sec_dis_num)=abs(P_amp_along_axial{sec_dis_num,section_plot}(3));
        
    end
    h1(section_plot)=plot(x_plot, p_plot,'Parent',hAxes(1),...
        'color','k',...
        'LineStyle','-',...
        'Linewidth',4,...
        'DisplayName',sprintf('Pressure'));
    hold on 
    if max(p_plot)>p_max
        p_max=max(p_plot);
    else
    end
    hold on 
    
end
% set fonts, labels, lims, ticks and so on
set(hAxes(1),'YColor','k','Box','on');
set(hAxes(1),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(1),'$x$ ($m$)','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(1),'$|\widetilde{p}|$ (Arbitrary scale)','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=min(x_all);
xlim_u=max(x_all);
xTick=xlim_d:(xlim_u-xlim_d)/5:xlim_u;
xTick_label=xTick;
set(hAxes(1),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.2f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=0;
ylim_u=p_max;
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(1),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.1f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(h1(1));
set(legend1,...
'Location','East',...
'visible','off',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
grid on
%********************
hAxes(2)=axes('Unit','pixels','position',[750 100 500 300]);
hold on
u_max=0;
for section_plot=1:1:CI.TP.numSection-1
    for sec_dis_num=1:1:length(x_section_discr)
        x_plot(sec_dis_num)=x_discr{sec_dis_num,section_plot}(1);
        p_plot(sec_dis_num)=abs(P_amp_along_axial{sec_dis_num,section_plot}(1));
        u_plot(sec_dis_num)=abs(P_amp_along_axial{sec_dis_num,section_plot}(3))/CI.TP.u_mean(2,section_plot);
        
    end
    h3(section_plot)=plot(x_plot, u_plot,'Parent',hAxes(2),...
        'color','k',...
        'LineStyle','--',...
        'Linewidth',4,...
        'DisplayName',sprintf('Velocity'));
    hold on 
    if max(u_plot)>u_max
        u_max=max(u_plot);
    else
    end
    
end
% set fonts, labels, lims, ticks and so on
set(hAxes(2),'YColor','k','Box','on');
set(hAxes(2),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(2),'$x$ ($m$)','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(2),'$|\widetilde{u}|/\bar{u}$ (Arbitrary scale)','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=min(x_all);
xlim_u=max(x_all);
xTick=xlim_d:(xlim_u-xlim_d)/5:xlim_u;
xTick_label=xTick;
set(hAxes(2),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.2f');
ylim_d=0;
ylim_u=u_max;
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(2),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.2f');
legend1=legend(h3(1));
set(legend1,...
'Location','East',...
'visible','off',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
grid on

%********************
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[1,1,1800,650])
set(gcf,'PaperPositionMode','auto');
%****************************
savename = ['Modeshape_test'];
% saveas(gcf,savename,'fig');
saveas(gcf,savename,'epsc2')
% print(savename,'-depsc')
savenameEps = [savename '.eps'];
savenamePdf = [savename '.pdf'];
eps2pdf(savenameEps,savenamePdf);
%-------------------------------end--------------------------------------
%
%----------------------Entropy convection transfer function ---------------
%
function Te = Fcn_TF_entropy_convection(s)
global CI
tau     = CI.setup.tau_Cs;
k       = CI.setup.ke;
Te = k*exp((tau*s)^2/4);
%
% ----------------------linear Flame transfer function --------------------
%         
function F = Fcn_flame_model(s,indexHA)
global CI
F = CI.setup.FM.a_f*exp(-CI.setup.FM.tau_f*s);     
%
% -----------------------------end-----------------------------------------