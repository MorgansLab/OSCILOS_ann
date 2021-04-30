% ----------------------------Calculate contour values---------------------
%
function [Freq_plot, GR_plot, F] = Fcn_contour_plot(Figure_num)
global CI
%
FreqMin=CI.EIG.Scan.FreqMin;
FreqMax=CI.EIG.Scan.FreqMax;
GRMin  =CI.EIG.Scan.GRMin;
GRMax  =CI.EIG.Scan.GRMax;
CI.EIG.Cont.FreqSp=linspace(FreqMin, FreqMax, 200); 
CI.EIG.Cont.GRSp  =linspace(GRMin, GRMax, 200); 
n   = length(CI.EIG.Cont.GRSp);
m   = length(CI.EIG.Cont.FreqSp);
F   = zeros(n,m);
for ss = 1:n
    GR = CI.EIG.Cont.GRSp(ss);
    for kk  = 1:m
        s           = GR+1i*2*pi*CI.EIG.Cont.FreqSp(kk);
        Freq_plot(ss,kk) = CI.EIG.Cont.FreqSp(kk);
        GR_plot(ss,kk)   = CI.EIG.Cont.GRSp(ss);                                                              
                F(ss,kk) = Fcn_DetEqn_Linear(s);                             
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot of the unstable modes in the s-plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure(Figure_num);
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(3).*(1/100) scrsz(4).*(1/20) scrsz(3)*5/7 scrsz(4).*(6/7)])
set(h,'name','Mode analysis for a choked-choked straight duct','numbertitle','off');
%************
hAxes(1)=axes('Unit','pixels','position',[150 100 500 500]);
hold on
contourf(GR_plot,Freq_plot,log10(abs(F)));
colormap hot;
hcb=colorbar;
plot(hAxes(1),CI.Eigenmode.GR ,CI.Eigenmode.Freq,'p',...
        'markersize',8,'color','k','markerfacecolor',[1,1,1])
set(hAxes(1),'YColor','k','Box','on');
set(hAxes(1),'FontName','Helvetica','FontSize',20,'LineWidth',1)
xlabel(hAxes(1),'Growth rate (1/s)','Color','k','Interpreter','LaTex','FontSize',20);
ylabel(hAxes(1),'Frequency (Hz)','Color','k','Interpreter','LaTex','FontSize',20);
set(hAxes(1),'ylim',[FreqMin FreqMax],'yTick',FreqMin:(FreqMax-FreqMin)/5:FreqMax,'YAxisLocation','left','Color','w');
set(hAxes(1),'xlim',[GRMin GRMax],'xTick',GRMin:(GRMax-GRMin)/4:GRMax);
grid on
%%%% add the grid on
hAxes(2)=axes('Unit','pixels','position',[100 100 500 500]);
hold on
set(hAxes(2),'ylim',[FreqMin FreqMax],'yTick',FreqMin:(FreqMax-FreqMin)/5:FreqMax,'yticklabel',[],...
    'YAxisLocation','left','Color','none');
set(hAxes(2),'xlim',[GRMin GRMax],'xTick',GRMin:(GRMax-GRMin)/4:GRMax,'xticklabel',[],...
    'xcolor','w','ycolor','w','gridlinestyle','-.');
set(hAxes(2),'FontName','Helvetica','FontSize',20,'LineWidth',1)
grid on
%% add colorbar
hAxes(3)=axes('Unit','pixels','position',[100 100 500 500]);
hold on;
set(hAxes(3),'ylim',[FreqMin FreqMax],'yTick',FreqMin:(FreqMax-FreqMin)/5:FreqMax,'yticklabel',[],...
    'YAxisLocation','left','Color','none');
set(hAxes(3),'xlim',[GRMin GRMax],'xTick',GRMin:(GRMax-GRMin)/4:GRMax,'xticklabel',[]);
set(hAxes(3),'FontName','H1elvetica','FontSize',20,'LineWidth',1)
set(hAxes(1),'Unit','pixels','position',[100 100 500 500]);
set(hcb,'Fontsize',20,'box','on','Unit','pixels')
set(hcb,'position',[620,100,25,500]);
%********************
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[1,1,1800,650])
set(gcf,'PaperPositionMode','auto');
%****************************
savename = ['Contour_plot'];
% saveas(gcf,savename,'fig');
saveas(gcf,savename,'epsc2')
% print(savename,'-depsc')
savenameEps = [savename '.eps'];
savenamePdf = [savename '.pdf'];
eps2pdf(savenameEps,savenamePdf);
% command = ['convert -interlace none -density 400 -quality 100 ' savename '.pdf ' savename '.png'];
% [status,cmdout] = system(command);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------------end--------------------------------------