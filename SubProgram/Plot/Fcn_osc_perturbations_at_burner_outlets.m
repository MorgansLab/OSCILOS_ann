function [A_u_d_complex, A_p_d_complex]=Fcn_osc_perturbations_at_burner_outlets(figure_num)

global CI
mode_solution=CI.EIG.nonlinear_solution_unique{1,1};
s=1i*mode_solution(1)+CI.fixed_growthrate; 
lambda(1,1:CI.setup.N)=mode_solution(2:(1+CI.setup.N))+1i*mode_solution((2+CI.setup.N):(1+2*CI.setup.N));
lambda(1,1+CI.setup.N)=mode_solution(2+2*CI.setup.N);
lambda(1,(2+CI.setup.N):(1+2*CI.setup.N))=mode_solution((3+2*CI.setup.N):(2+3*CI.setup.N))+1i*mode_solution((3+3*CI.setup.N):(2+4*CI.setup.N));
ss=find(CI.setup.DuctIndex==1);
ss_premixer=ss(end);

%% Perturbations just before the flames
[p_bf_n,u_bf_n]=Fcn_p_u_bf(s,ss_premixer);  % To calculate oscillations before the flames.
p_bf_n=p_bf_n.*lambda(1,:);  
u_bf_n=u_bf_n.*lambda(1,:);            % This is where lambda comes in

%% Heat release perturbations
N_all=2*CI.setup.N+1;

indexHA=1;                              % Only one flame interface is considered in the present version of OSCILOS_Annular
D=CI.setup.Premixer_Num;
Q_mean_d=CI.TP.Q(indexHA)/D;            % This mean heat release rate is assumed to be the same for all premixers.
[p_3_n,u_3_n]=Fcn_u_3_n(s,ss_premixer);  % To calculate oscillations before the flames.
p_3_n=p_3_n.*lambda(1,:);  
u_3_n=u_3_n.*lambda(1,:);            % This is where lambda comes in
u_mean_d=CI.TP.u_mean(1,ss_premixer);   % This mean velocity is assumed to be the same for all premixers. This is inside the premixers 
n_all=-CI.setup.N:1:CI.setup.N;
D=CI.setup.Premixer_Num;
phi_d_all_coe=0:1/D:(D-1)/D;
phi_d_all=2*pi*phi_d_all_coe;
for d=1:1:D
    A_p_d_complex(d,1) = sum(p_3_n.*exp(1i.*n_all.*phi_d_all(d)));
    A_u_d_complex(d,1) = sum(u_3_n.*exp(1i.*n_all.*phi_d_all(d)))/u_mean_d;
end
A_u_d(:,1)=abs(A_u_d_complex(:,1));
A_p_d(:,1)=abs(A_p_d_complex(:,1));

for n_scan=1:1:N_all
    n0=n_all(n_scan);

for d=1:1:D
    A_u_d_Q(d) = sum(u_3_n.*exp(1i.*n_all.*phi_d_all(d)))/u_mean_d;
    Flame_model_d(d)=Fcn_flame_model(s,A_u_d_Q(d),indexHA);
    % Q_d
    Q_d(d)   = Q_mean_d*Flame_model_d(d)*A_u_d_Q(d);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Modeshape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure(figure_num);
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(3).*(1/100) scrsz(4).*(1/20) scrsz(3)*6.5/7 scrsz(4).*(6/7)])
set(h,'name','Mode analysis for a choked-choked straight duct','numbertitle','off');
%************
hAxes(1)=axes('Unit','pixels','position',[150 100 700 250]);
hold on
hA1=bar(A_u_d(:,1));
% set fonts, labels, lims, ticks and so on
set(hAxes(1),'YColor','k','Box','on');
set(hAxes(1),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(1),'$d$','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(1),'$|\widetilde{u}_d|/\bar{u}$','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=0;
xlim_u=17;
xTick=1:1:16;
xTick_label=xTick;
set(hAxes(1),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.0f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=0;%min(A_u_d(:,3+1));
ylim_u=max(A_u_d(:,1));
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(1),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.2f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(hA1(1));
set(legend1,...
'Location','East',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
legend(hAxes(1),'off');
grid on
%************
hAxes(4)=axes('Unit','pixels','position',[150 430 700 250]);
hold on
hA4=bar(A_p_d(:,1));
% set fonts, labels, lims, ticks and so on
set(hAxes(4),'YColor','k','Box','on');
set(hAxes(4),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(4),'$d$','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(4),'$|\widetilde{p}|~(Pa)$','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=0;
xlim_u=17;
xTick=1:1:16;
xTick_label=xTick;
set(hAxes(4),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.0f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=0;%min(A_p_d(:,3+1));
ylim_u=max(A_p_d(:,1));
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(4),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.2f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(hA4(1));
set(legend1,...
'Location','East',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
legend(hAxes(4),'off');
grid on
%%%%%%%%%%%%%%%%%%%%%%%%
hAxes(5)=axes('Unit','pixels','position',[1000 100 700 250]);
hold on
hA1=bar(abs(Q_d));
% set fonts, labels, lims, ticks and so on
set(hAxes(5),'YColor','k','Box','on');
set(hAxes(5),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(5),'$d$','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(5),'$|\widetilde{Q}_d|~(J/s)$','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=0;
xlim_u=17;
xTick=1:1:16;
xTick_label=xTick;
set(hAxes(5),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.0f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=0;%min(A_u_d(:,3+1));
ylim_u=max(abs(Q_d));
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(5),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.2f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(hA1(1));
set(legend1,...
'Location','East',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
legend(hAxes(5),'off');
grid on
%********************
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[1,1,1800,650])
set(gcf,'PaperPositionMode','auto');
%****************************
savename = ['Perturbations_at_burner_outlet'];
% saveas(gcf,savename,'fig');
saveas(gcf,savename,'epsc2')
% print(savename,'-depsc')
savenameEps = [savename '.eps'];
savenamePdf = [savename '.pdf'];
eps2pdf(savenameEps,savenamePdf)