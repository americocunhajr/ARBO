

% plot MC convergence metric
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' number of MC realizations';
ylab   = ' convergence metric';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
%ymin   = 0.22;
%ymax   = 0.25;
gname  = [num2str(case_name),'__MC_conv'];
flag   = 'eps';
fig1   = graph_type1x((1:Ns),MC_conv,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% -----------------------------------------------------------



samp_pts = randi([1 Ns],1,5);


% plot cumulative confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' nominal';
leg5   = ' data';
xlab   = ' time (weeks)';
ylab   = ' number of people';
xmin   = t0/7;
xmax   = t1/7;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_C'];
gname1 = [num2str(case_name),'__stat_C'];
flag   = 'eps';
% fig2a  = graph_ci1N(time,MC_Qdisp(samp_pts,:),...
%                     Qdisp_upp,...
%                     Qdisp_low,...
%                     gtitle,leg1,leg2,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

% fig2a  = graph_ci2(time/7,C_smp_avg,...
%                    time/7,C_std,...
%                    C_upp,...
%                    C_low,...
%                    gtitle,leg3,leg4,leg2,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

fig2a  = graph_ci2data(time/7,C_smp_avg,...
                       time/7,Nominal_C,...
                       C_upp,C_low,...
                       day/7,CData,...
                       gtitle,leg3,leg4,leg2,leg5,...
                       xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

% fig2aa = graph_type2(time,C_smp_avg,...
%                      time,C_std,...
%                      gtitle,leg3,leg4,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);

%close(fig2a);
%close(fig2aa);
% -----------------------------------------------------------


% plot new cases confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' nominal';
leg5   = ' data';
xlab   = ' time (weeks)';
ylab   = ' number of people';
xmin   = 1;
xmax   = Nweeks;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_NewCases'];
gname1 = [num2str(case_name),'__stat_NewCases'];
flag   = 'eps';
% fig2b  = graph_ci1N(time,MC_Qdisp(samp_pts,:),...
%                     Qdisp_upp,...
%                     Qdisp_low,...
%                     gtitle,leg1,leg2,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

%fig2b  = graph_ci2((1:Nweeks),NewCases_smp_avg,...
%                   (1:Nweeks),NewCases_std,...
%                   NewCases_upp,...
%                   NewCases_low,...
%                   gtitle,leg3,leg4,leg2,...
%                   xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

fig2b  = graph_ci2data((1:Nweeks),NewCases_smp_avg,...
                       (1:Nweeks),Nominal_NewCases,...
                       NewCases_upp,NewCases_low,...
                       (1:Nweeks),NewCasesData,...
                       gtitle,leg3,leg4,leg2,leg5,...
                       xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);


% fig2bb = graph_type2((1:Nweeks),NewCases_smp_avg,...
%                      (1:Nweeks),NewCases_std,...
%                      gtitle,leg3,leg4,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);

%close(fig2b);
%close(fig2bb);
% -----------------------------------------------------------



% % plot displacement time-averaged PDF
% % -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' time averaged displacement';
% ylab   = ' probability density function';
% %xmin   = -5.0;
% %xmax   =  5.0;
% %ymin   = 0.0;
% %ymax   = 2.5;
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__pdf_disp'];
% flag   = 'eps';
% fig3a  = graph_bar_curve1(Qdisp_bins,Qdisp_freq,...
%                           Qdisp_supp,Qdisp_ksd,gtitle,...
%                           xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
% %close(fig3a);
% % -----------------------------------------------------------

% 
% % plot power PDF at time t
% % -----------------------------------------------------------
% tt1 = round(0.25*Ndt);
% tt2 = round(0.50*Ndt);
% tt3 = round(0.75*Ndt);
% tt4 = Ndt;
% 
% %gtitle = ' tower horizontal displacement';
% gtitle = ' ';
% xlab   = ' power';
% ylab   = ' probability density function';
% xmin   = 0.0;
% xmax   = 0.06;
% ymin   = 0.0;
% ymax   = 200.0;
% %xmin   = 'auto';
% %xmax   = 'auto';
% %ymin   = 'auto';
% %ymax   = 'auto';
% gname1  = [num2str(case_name),'__pdf_power_',num2str(tt1)];
% gname2  = [num2str(case_name),'__pdf_power_',num2str(tt2)];
% gname3  = [num2str(case_name),'__pdf_power_',num2str(tt3)];
% gname4  = [num2str(case_name),'__pdf_power_',num2str(tt4)];
% flag   = 'eps';
% fig4a1 = graph_bar_curve1(power_bins_t(:,tt1),power_freq_t(:,tt1),...
%                           power_supp_t(:,tt1),power_ksd_t(:,tt1),gtitle,...
%                           xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
% fig4a2 = graph_bar_curve1(power_bins_t(:,tt2),power_freq_t(:,tt2),...
%                           power_supp_t(:,tt2),power_ksd_t(:,tt2),gtitle,...
%                           xlab,ylab,xmin,xmax,ymin,ymax,gname2,flag);
% fig4a3 = graph_bar_curve1(power_bins_t(:,tt3),power_freq_t(:,tt3),...
%                           power_supp_t(:,tt3),power_ksd_t(:,tt3),gtitle,...
%                           xlab,ylab,xmin,xmax,ymin,ymax,gname3,flag);
% fig4a4 = graph_bar_curve1(power_bins_t(:,tt4),power_freq_t(:,tt4),...
%                           power_supp_t(:,tt4),power_ksd_t(:,tt4),gtitle,...
%                           xlab,ylab,xmin,xmax,ymin,ymax,gname4,flag);
% %close(fig4a1);
% %close(fig4b1);
% %close(fig4c1);
% %close(fig4d1);
% % -----------------------------------------------------------
