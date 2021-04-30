function [] = Fcn_nonlinear_mode_evolution(S_nonlinear, Figure_num)
% This is to plot the nonlinear evolution of a chosen mode
global CI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results     =S_nonlinear;
Growth_rate =Results(:,1);
Frequency   =Results(:,2)/2/pi;
Mode_num    =-CI.setup.N:1:CI.setup.N;

for nn=1:1:(CI.setup.N)
    lambda(:,nn)=Results(:,2+nn) + 1i*Results(:,2+CI.setup.N+nn); 
    lambda(:,nn+CI.setup.N+1)=Results(:,3+2*CI.setup.N+nn) + 1i*Results(:,3+3*CI.setup.N+nn);
end
lambda(:,CI.setup.N+1)=Results(:,3+2*(CI.setup.N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Modeshape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Line setting parameters
Number=CI.setup.N;


delta_color=5/(12*Number);
color(1:(2*Number+1),1)=[0:delta_color:delta_color*2*Number].';
color(1:(2*Number+1),2)=[0:delta_color:delta_color*2*Number].';
color(1:(2*Number+1),3)=[0:delta_color:delta_color*2*Number].';
delta_width=2;
width(1:(2*Number+1))=delta_width*[1:1:(2*Number+1)];



h=figure(Figure_num);
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(3).*(1/100) scrsz(4).*(1/20) scrsz(3)*5/7 scrsz(4).*(6/7)])
set(h,'name','Mode analysis for a choked-choked straight duct','numbertitle','off');
%************
hAxes(1)=axes('Unit','pixels','position',[150 100 700 200]);
hold on
h1=plot(Growth_rate,Frequency,...
        'color','k',...
        'LineStyle','--',...
        'Linewidth',4,...
        'DisplayName',sprintf('Frequency'));
% set fonts, labels, lims, ticks and so on
set(hAxes(1),'YColor','k','Box','on');
set(hAxes(1),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(1),'Growth rate ($s^{-1}$)','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(1),'Frequency ($Hz$)','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=min(Growth_rate);
xlim_u=max(Growth_rate);
xTick=xlim_d:(xlim_u-xlim_d)/5:xlim_u;
xTick_label=xTick;
set(hAxes(1),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.0f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=min(Frequency);
ylim_u=max(Frequency);
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(1),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.0f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(h1(1));
set(legend1,...
'Location','East',...
'FontSize',20,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
legend(hAxes(1), 'off')
grid on
%********************
hAxes(2)=axes('Unit','pixels','position',[150 380 700 350]);
hold on
Modes_plot=2:1:(2*CI.setup.N+1-1);

for nn=1:1:length(Modes_plot)
    h2(nn)=plot(Growth_rate,abs(lambda(:,Modes_plot(nn))),...
        'color',color(Modes_plot(nn),:),...
        'LineStyle','--',...
        'Linewidth',width(Modes_plot(nn)),...
        'DisplayName',sprintf('n=%d', Mode_num(Modes_plot(nn))));
end
% set fonts, labels, lims, ticks and so on
set(hAxes(2),'YColor','k','Box','on');
set(hAxes(2),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(2),'Growth rate ($s^{-1}$)','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(2),'Mode strength ($Pa$)','Color','k','Interpreter','LaTex','FontSize',20);
xlim_d=min(Growth_rate);
xlim_u=max(Growth_rate);
xTick=xlim_d:(xlim_u-xlim_d)/5:xlim_u;
xTick_label=xTick;
set(hAxes(2),'xlim',[xlim_d xlim_u],'xTick',xTick);
xtickformat('%.0f');
%set(gca,'XTickLabel',sprintf('%1.2f|',xTick_label), 'Fontsize',20);
ylim_d=min(min(abs(lambda(:,:))));
ylim_u=max(max(abs(lambda(:,:))));
yTick=ylim_d:(ylim_u-ylim_d)/5:ylim_u;
set(hAxes(2),'ylim',[ylim_d ylim_u],'yTick',yTick);
ytickformat('%.2f');
%set(gca,'YTickLabel',sprintf('%1.2f',yTick), 'Fontsize',20);
legend1=legend(h2);
set(legend1,...
'Location','best',...
'FontSize',18,...
'Interpreter','LaTex');
set(legend1, 'Box', 'off');
grid on
%********************
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[1,1,1800,650])
set(gcf,'PaperPositionMode','auto');
%****************************
savename = ['Mode_evolution_test'];
% saveas(gcf,savename,'fig');
saveas(gcf,savename,'epsc2')
% print(savename,'-depsc')
savenameEps = [savename '.eps'];
savenamePdf = [savename '.pdf'];
eps2pdf(savenameEps,savenamePdf)

% command = ['convert -interlace none -density 300 -quality 80 ' savename '.pdf ' savename '.png'];
% [status,cmdout] = system(command);

% % -------------------------------end--------------------------------------