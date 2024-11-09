%% 
% x0 = [1;-0.375;4];
% xf = [-0.1,0.1,-1];
% x0 = [0;0;3];
% x0 = [-0.5;1;-3];
% xf = [0;0;0];
x0 = [1;-0.375;3.999];
xf = [0;0;4];
%% 
M_max = [1;1;1.5;4];
M_min = [-1;-1;-1.5;-4];
tic
[orders,signs,tangents,arctimes] = plan_nth_order(x0,xf,M_max,M_min,true,0,1e-6);
toc
%% 
Ts = 1e-3;
[xs,ts] = interpolate_MIM(x0,orders,signs,tangents,arctimes,M_max(1),M_min(1),Ts,0,true);
%% 
ylabels = ["Jerk","Acceleration","Velocity","Position"];
figure
for i = 1:4
    subplot(4,1,i)
    plot(ts,xs(i,:),'LineWidth',2)
    ax = gca;
    ax.FontName = "Times New Roman";
%     ax.FontSize = 16;
    ax.LineWidth = 1;
    ylabel(ax,ylabels(i))
    ax.XGrid = "on";
    ax.YGrid = "on";
    ax.XLim = ts([1,end]);
end
xlabel(ax,"Time")