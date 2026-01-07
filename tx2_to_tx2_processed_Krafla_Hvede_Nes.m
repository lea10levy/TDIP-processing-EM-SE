%% Krafla
% addpath('Data_Krafla\')
% addpath('mat\')
% path_out='Data_Krafla\';
% 
% name_SE_fit_Krafla_EM='Krafla_SEfit_log_curvatureMax_Min_End_deg5_v2.mat';
% load(name_SE_fit_Krafla_EM,'num_depth_ref','remaining_gates','remaining_m','xdata0_excl_EM','ydata0_excl_EM','val_clog_peak_all','all_locs_all','time_curvature_max_all','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2','coeff_alpha','coeff_tauKWW', 'coeff_beta','Mapp_cell','IP_flag_cell','t_center','current_cell','my_voltage','Rho_save','Res_save','file_in')
% manual_or_not=1;
% transpose_or_not=0;

%% Hvedermaken
addpath('Data_Hvedemarken\')
addpath('mat\')
path_out='Data_Hvedemarken\';
name_SE_fit_Hvede_EM='Hvede_SEfit_log_curvatureMax_Min_End_deg5.mat';
load(name_SE_fit_Hvede_EM,'num_depth_ref','remaining_gates','remaining_m','xdata0_excl_EM','ydata0_excl_EM','val_clog_peak_all','all_locs_all','time_curvature_max_all','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2','coeff_alpha','coeff_tauKWW', 'coeff_beta','Mapp_cell','IP_flag_cell','t_center','current_cell','my_voltage','Rho_save','Res_save','file_in')
manual_or_not=1;
transpose_or_not=0;

%% Nesjavellir
% addpath('Data_Nesjavellir\')
% addpath('mat\')
% path_out='Data_Nesjavellir\';
% name_SE_fit_Nes_EM='Nes_SEfit_log_curvatureMax_Min_End_deg5_v3.mat';
% load(name_SE_fit_Nes_EM,'my_proc_time','num_depth_ref','remaining_gates','remaining_m','xdata0_excl_EM','ydata0_excl_EM','val_clog_peak_all','all_locs_all','time_curvature_max_all','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2','coeff_alpha','coeff_tauKWW', 'coeff_beta','Mapp_cell','IP_flag_cell','t_center','current_cell','my_voltage','Rho_save','Res_save','file_in')
% manual_or_not=0;
% transpose_or_not=1;

%% load comparison data
% name_comp='Krafla_manual_auto_deg5_normr1_v2.mat';
% load(name_comp,'comp_dp','count_both_kept','count_both_removed','count_discr','auto_kept_not_manual','count_auto_kept_not_manual','manual_kept_notauto','comp_gates','count_gates_both_kept','count_gates_both_removed','count_gates_discr')

name_comp='Hvede_manual_auto_deg5_normr1_v2.mat';
load(name_comp,'comp_dp','count_both_kept','count_both_removed','count_discr','auto_kept_not_manual','count_auto_kept_not_manual','manual_kept_notauto','comp_gates','count_gates_both_kept','count_gates_both_removed','count_gates_discr')


%% find the flagged gates in automatic processing
num_dataset=numel(file_in);
% num_depth_ref =1678;
id_not_flagged_EM=cell(num_dataset,num_depth_ref);
is_dp_flag_manual=NaN(num_dataset,num_depth_ref);
is_dp_flag_auto=NaN(num_dataset,num_depth_ref);
% is_gate_flag_manual=cell(num_dataset,num_depth_ref);
% is_gate_flag_auto=cell(num_dataset,num_depth_ref);
% id_not_flagged_noEM=cell(num_dataset,num_depth_ref);
bloc_flagged_all=cell(num_dataset,1);
num_dp=NaN(num_dataset,1);
NGates=NaN(num_dataset,1);
removed_after_step_7=NaN(num_dataset,num_depth_ref);
removed_after_step_6=NaN(num_dataset,num_depth_ref);
removed_before_step_6=NaN(num_dataset,num_depth_ref);
% res_thres=0.2;
do_plot=1;
% rel_res_thres=0.1;
for k=1:num_dataset
    num_dp(k)=size(Mapp_cell{k},1);
    NGates(k)=size(Mapp_cell{k},2);
    bloc_flags=zeros(num_dp(k),NGates(k));
    id_flagged=ones(1,NGates(k)); 
    for i=1:num_dp(k)
        if isempty(gof_se_2{k,i})
            fitresult_se_EM=fitresult_se{k,i};
            gof_se_EM=gof_se{k,i};
            output_se_EM=output_se{k,i};
            se_fit_mapp_EM=se_fit_mapp{k,i};
            xdata_remaining=remaining_gates{k,i};
            ydata_remaining=remaining_m{k,i};
%             fprintf('first \n')
        else
            fitresult_se_EM=fitresult_se_2{k,i};
            gof_se_EM=gof_se_2{k,i};
            output_se_EM=output_se_2{k,i};
            se_fit_mapp_EM=se_fit_mapp_2{k,i};
            xdata_remaining=xdata0_excl_EM{k,i};
            ydata_remaining=ydata0_excl_EM{k,i};
%             fprintf('second \n')
        end
        xdata0_tmp = t_center{k}(1:end); 
        ydata0_tmp = Mapp_cell{k}(i,1:end)';

       % flags in manual proc
       if manual_or_not
            ind_gate_non_flagged_original=find(IP_flag_cell{k}(i,:)<1);
            xdata0_noflag = xdata0_tmp(ind_gate_non_flagged_original);
            ydata0_noflag = ydata0_tmp(ind_gate_non_flagged_original);
             if isempty(ydata0_noflag)
                is_dp_flag_manual(k,i)=1;
            else
                is_dp_flag_manual(k,i)=0;
            end
       end
        
        if isempty(fitresult_se_EM)
            warning(['no SE fit for DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=1;
            removed_before_step_6(k,i)=1;
        elseif gof_se_EM.rsquare<0.99
            warning(['Bad SE fit for DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=1;
            removed_before_step_6(k,i)=0;
        elseif coeff_alpha(k,i)>1000
            warning(['alpha too high for DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=0;
            removed_before_step_6(k,i)=0;
        elseif coeff_beta(k,i)<0.2
            warning(['beta too low for DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=0;
            removed_before_step_6(k,i)=0;
        elseif coeff_beta(k,i)==1
            warning(['beta too highfor DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=0;
            removed_before_step_6(k,i)=0;
        elseif coeff_tauKWW(k,i)<1
            warning(['Tau too low for DP =' num2str(i) 'm for file' file_in{k}])
            id_not_flagged_EM{k,i}=zeros(1,NGates(k)); % 1 if not flagged
            EM_ok=0;
            is_dp_flag_auto(k,i)=1;
            removed_after_step_7(k,i)=1;
            removed_after_step_6(k,i)=0;
            removed_before_step_6(k,i)=0;
        else
            EM_ok=1;      
            % add the transpose for Nes, but not for Hvede/Krafla
            if transpose_or_not
                id_not_flagged_EM{k,i}=ismember(xdata0_tmp',xdata_remaining); % 1 if not flagged
            else
                id_not_flagged_EM{k,i}=ismember(xdata0_tmp,xdata_remaining); % 1 if not flagged
            end
            is_dp_flag_auto(k,i)=0;
            removed_after_step_7(k,i)=0;
            removed_after_step_6(k,i)=0;
            removed_before_step_6(k,i)=0;
        end
        
        % make table comparison auto vs manual proc
        
        % plot section
    
%         if ismember(k,[5,6,7]) %ISL10-11-12
        if (do_plot) 
%             if manual_kept_notauto(k,i)
%             if auto_kept_not_manual(k,i)
            if comp_dp(k,i)==0
                clf;
%             if comp_dp(k,i)==0
%         if is_dp_flag_manual(k,i)<is_dp_flag_auto(k,i) %auto_kept_not_manual(k,i)%manual_kept_notauto(k,i)    
%             if ismember(i,[4,7,10,11,12,16,117,173,728,1383,1396,1397,1402,1408,1410,1413,1417,1418,1422,1433,1445,1495,1502,1571,1656,1669,1673,1675])
                lgd=[];
                n_leg=1;
%                 plot(xdata0_tmp,ydata0_tmp,'o','Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',9);
                loglog(xdata0_tmp, abs(ydata0_tmp),'or','MarkerFaceColor','r','MarkerSize',6)
                hold on
                loglog(xdata0_tmp, ydata0_tmp,'ob','MarkerFaceColor','b','MarkerSize',6)
                lgd{n_leg}='Unprocessed data - negative';
                n_leg=n_leg+1;
                lgd{n_leg}='Unprocessed data - positive';
                hold on
%                 if EM_ok
                    n_leg=n_leg+1;
                    plot(xdata_remaining,ydata_remaining,'o','Color','k','MarkerFaceColor','k','MarkerSize',8);
                    lgd{n_leg}='Data after EM+SE AutoProc';
                    if isempty(se_fit_mapp_EM)
                        txt_EM='No SE fit';
                    else
                        plot(xdata_remaining,se_fit_mapp_EM,'--','Color','k','LineWidth',2);
                        n_leg=n_leg+1;
                        lgd{n_leg}='SE fit';
                        txt_EM=['R^2_{SE} =' num2str(gof_se_EM.rsquare,4)];
                        if num2str(gof_se_EM.rsquare)>0.99
                            txt_SE_param=['\alpha =' num2str(round(coeff_alpha(k,i),2)) 'mV/V; \beta =' num2str(coeff_beta(k,i),1) '; \tau_{KWW} =' num2str(round(coeff_tauKWW(k,i),2)) ' ms'];
                            hText3=text(20,4000,txt_SE_param);
                            set(hText3,'FontName','Times New Roman','FontSize',12);
                        end
                    end
                    hText2=text(500,1000,txt_EM);
                    
%                     plot(remaining_gates{k,i},se_fit_mapp{k,i},'--','Color','b','LineWidth',2);
                    set(hText2,'FontName','Times New Roman','FontSize',12);
%                     n_leg=n_leg+1;
%                     lgd{n_leg}='SE fit 1';

%                 end
                if is_dp_flag_manual(k,i)==0
                    plot(xdata0_noflag,ydata0_noflag,'*','Color',[0 1 0],'MarkerSize',3,'LineWidth',1);
                    n_leg=n_leg+1;
                    lgd{n_leg}='Manual processing';
                    plot(xdata0_noflag,abs(ydata0_noflag),'*','Color',[0 1 0],'MarkerSize',3,'LineWidth',1);
                end
                
%                 txt_title=[file_in{k}(1:end-4) '- DP : ' num2str(i)]; % Krafla
                txt_title=[file_in{k}(1:12) '- DP : ' num2str(i)]; % Hvede
                hTitle=title(txt_title,'Interpreter','none');
                
                hLegend=legend(lgd,'Location','SouthWest');
                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
        
                hXLabel=xlabel('Time (msec)');
                hYLabel=ylabel('Polarizability \eta (mV/V)');
                xlim([1e0,1e4])
                ylim([1e-2,1e4])
                grid on;
                hold off
                set(gca,'FontName','Times New Roman','FontSize',14);
                set([hTitle,hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
%                 name_out=['png\Krafla_both_kept_deg5_normr1_v2\' file_in{k}(1:end-4) '_DP' num2str(i) '_auto_EM_manual.png'];

%                 name_out=['png\Hvede_manual_kept_not_auto_deg5_normr1_v2\' file_in{k}(1:end-4) '_DP' num2str(i) '_auto_EM_manual.png'];
%                 name_out=['png\Krafla_auto_kept_not_manual_deg5_normr1_v2\' file_in{k}(1:end-4) '_DP' num2str(i) '_auto_EM_manual.png'];
                name_out=['pdf\Hvede_bothkept\' file_in{k}(1:12) '_DP' num2str(i) '_auto_EM_manual.pdf'];
                exportgraphics(gcf,name_out,'Resolution',300)
                name_out=['png\Hvede_bothkept\' file_in{k}(1:12) '_DP' num2str(i) '_auto_EM_manual.png'];
                exportgraphics(gcf,name_out,'Resolution',300)
            end
        end
        if i<num_dp(k)
            bloc_flags(i,:)=id_flagged-id_not_flagged_EM{k,i};
        end
    end
    bloc_flagged_all{k}=bloc_flags;
end

%% count kept / removed at different steps
count_kept_after_step_7=sum(removed_after_step_7==0,'all');
count_kept_after_step_6=sum(removed_after_step_6==0,'all');
count_kept_before_step_6=sum(removed_before_step_6==0,'all');


%% make table comparison auto vs manual
% attention tout était inversé sur both kept/both removed (dp et gates),
% maintenant j'ai changé

comp_dp=is_dp_flag_auto+is_dp_flag_manual;
count_both_kept=sum(comp_dp==0,'all');
count_both_removed=sum(comp_dp==2,'all');
count_discr=sum(comp_dp==1,'all');
auto_kept_not_manual=is_dp_flag_auto<is_dp_flag_manual;
count_auto_kept_not_manual=sum(auto_kept_not_manual==1,'all');

manual_kept_notauto=is_dp_flag_auto>is_dp_flag_manual;
count_manual_kept_notauto=sum(manual_kept_notauto==1,'all');

comp_gates=cell(num_dataset,1);
count_gates_both_kept=cell(num_dataset,1);
count_gates_both_removed=cell(num_dataset,1);
count_gates_discr=cell(num_dataset,1);
for k=1:num_dataset
    comp_gates{k}=bloc_flagged_all{k}+IP_flag_cell{k};
    count_gates_both_kept{k}=sum(comp_gates{k}==0,'all'); % il faudrait plutot faire un histogramme par datapoint ici
    count_gates_both_removed{k}=sum(comp_gates{k}==2,'all');
    count_gates_discr{k}=sum(comp_gates{k}==1,'all');
end

name_final='Krafla_manual_auto_deg5_normr1_v2.mat';
save(name_final,'comp_dp','count_both_kept','count_both_removed','count_discr','auto_kept_not_manual','count_auto_kept_not_manual','manual_kept_notauto','comp_gates','count_gates_both_kept','count_gates_both_removed','count_gates_discr')


%% read old tx2 and replace the flags
for k=1:n
    fileOut=[path_out file_in{k}(1:end-4) '_automatic.tx2']; 

    fid=fopen(file_in{k},'r');  headerline = fgetl(fid);fclose(fid);
    tx2_table=dlmread(file_in{k},'\t',1,0);
    num_dp_tmp=numel(tx2_table(:,1));

    num_gates=tx2_table(1,25);
    id_mapp=26:26+num_gates-1;
    id_mdly=26+num_gates;
    id_gate_width=26+num_gates+1:26+num_gates*2;
    id_mflag=26+num_gates*3+1:26+num_gates*4;
    id_current=26+num_gates*4+2;
    Mapp=tx2_table(:,id_mapp);
    Gate_width=tx2_table(:,id_gate_width);

    tx2_table(:,id_mflag)=bloc_flagged_all{k};

    I_empty = [];
    headerwords = split(headerline);
    for i=1:numel(headerwords)
        if isempty(headerwords{i})
            I_empty = [I_empty i];
        end
    end
    headerwords(I_empty) = [];
    fid=fopen(fileOut,'w');
    for i = 1:numel(headerwords)
        fprintf(fid,'%10s\t',headerwords{i});
    end
    fprintf(fid,'\n');
    fclose(fid);
    dlmwrite(fileOut,tx2_table,'-append','delimiter','\t','precision','%10.5f');
end