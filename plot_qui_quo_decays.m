addpath('Data_Krafla\')

delimiterIn = ' ';
headerlinesIn = 1;

file_gates='ISL10_gates.txt';
gate_file=importdata(file_gates,delimiterIn,headerlinesIn);
data_gates=gate_file.data;
NGates=size(data_gates,1);
Gate_Center=data_gates(:,3);

file_qui='ISL10_MPA_qui_data_table.txt';
data_file_qui=importdata(file_qui,delimiterIn,headerlinesIn);
data_qui=data_file_qui.data;
col_headers_qui=data_file_qui.colheaders;
Ndata_qui=size(data_qui,1);
Ncol_qui=size(data_qui,2);
el_ABMN_qui=[data_qui(:,5) data_qui(:,7) data_qui(:,10) data_qui(:,11)];
flagIP_qui = data_qui(:,60:97);
IP_Mapp_qui = data_qui(:,104:141);
% IP_Mapp_qui_notflagged=IP_Mapp_qui(flagIP_qui==0);

file_quo='ISL10_MPA_quo_data_table.txt';
data_file_quo=importdata(file_quo,delimiterIn,headerlinesIn);
data_quo=data_file_quo.data;
col_headers_quo=data_file_quo.colheaders;
Ndata_quo=size(data_quo,1);
Ncol_quo=size(data_quo,2);
el_ABMN_quo=[data_quo(:,5) data_quo(:,7) data_quo(:,10) data_quo(:,11)];
flagIP_quo = data_quo(:,60:97);
IP_Mapp_quo = data_quo(:,104:141);
% IP_Mapp_quo_notflagged=IP_Mapp_quo(flagIP_quo==0);

num_dp=Ndata_qui;
clf;
% my_ABMN1=[13,27,20,21];
% my_ABMN2=[14,28,20,21];
% my_ABMN2=[1,29,19,21];
% my_ABMN1=[7,31,15,23];
% my_ABMN2=[8,22,15,16];
% my_ABMN1=[9,23,15,16];
my_ABMN2=[7,21,15,16];
my_ABMN1=[10,24,15,16];
[Lia1,Locb1]=ismember(el_ABMN_qui,my_ABMN1,'rows');
[Lia2,Locb2]=ismember(el_ABMN_qui,my_ABMN2,'rows');
idx=1:num_dp;
my_loc1=idx(Lia1);
my_loc2=idx(Lia2);
for i=[my_loc1,my_loc2]%:num_dp
    
    lgd=[];
    n_leg=1;
    
%     loglog(Gate_Center*1000, abs(IP_Mapp_qui(i,:)),'o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',4)
    lgd{n_leg}='Data excluded';
    hold on
    loglog(Gate_Center(flagIP_qui(i,:)==0)*1000, abs(IP_Mapp_qui(i,flagIP_qui(i,:)==0)),'or','MarkerFaceColor','r','MarkerSize',6)
    n_leg=n_leg+1;
    lgd{n_leg}='Data - negative';
%     hold on
    loglog(Gate_Center(flagIP_qui(i,:)==0)*1000, IP_Mapp_qui(i,flagIP_qui(i,:)==0),'ob','MarkerFaceColor','b','MarkerSize',6)    
    n_leg=n_leg+1;
    lgd{n_leg}='Data - positive';
    
    loglog(Gate_Center*1000, abs(IP_Mapp_quo(i,:)),'--k','LineWidth',2)
    n_leg=n_leg+1;
    lgd{n_leg}='Forward - negative';
%     hold on
    loglog(Gate_Center*1000, IP_Mapp_quo(i,:),'-k','LineWidth',2)    
    n_leg=n_leg+1;
    lgd{n_leg}='Forward - positive';

%     txt_title=[file_qui(1:end-4) '- DP : ' num2str(i)]; % Krafla
%     hTitle=title(txt_title,'Interpreter','none');
% 
% %     hLegend=legend(lgd,'Location','BestOutside');
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
% 
%     hXLabel=xlabel('Time (msec)');
%     hYLabel=ylabel('Polarizability \eta (mV/V)');
%     xlim([1e0,1e4])
%     ylim([1e-2,1e4])
%     grid on;
% %     hold off
%     set(gca,'FontName','Times New Roman','FontSize',14);
% %     set([hTitle,hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
%     set([hTitle,hXLabel,hYLabel],'FontName','Times New Roman','FontSize',14);
%     name_out=['png\Krafla_qui_quo\' file_qui(1:end-4) '_DP' num2str(i) '_qui_quo.png'];
%     exportgraphics(gcf,name_out,'Resolution',300)
%     name_out=['png\Krafla_qui_quo\' file_qui(1:end-4) '_DP' num2str(i) '_qui_quo.pdf'];
%     exportgraphics(gcf,name_out,'Resolution',300)
end

txt_title=[file_qui(1:end-4) '- El : ' num2str([my_ABMN1,my_ABMN2])]; % Krafla
hTitle=title(txt_title,'Interpreter','none');

%     hLegend=legend(lgd,'Location','BestOutside');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hXLabel=xlabel('Time (msec)');
hYLabel=ylabel('Polarizability \eta (mV/V)');
xlim([1e0,1e4])
ylim([1e-2,1e4])
grid on;
    hold off
set(gca,'FontName','Times New Roman','FontSize',14);
%     set([hTitle,hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
set([hTitle,hXLabel,hYLabel],'FontName','Times New Roman','FontSize',14);
name_out=['png\Krafla_qui_quo\' file_qui(1:end-4) '_Overlap_qui_quo_nogrey_3.png'];
exportgraphics(gcf,name_out,'Resolution',300)
name_out=['png\Krafla_qui_quo\' file_qui(1:end-4) '_Overlap_qui_quo_nogrey_3.pdf'];
exportgraphics(gcf,name_out,'Resolution',300)
