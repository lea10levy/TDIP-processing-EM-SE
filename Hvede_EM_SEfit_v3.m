addpath('Data_Hvedemarken\')
path_out='Data_Hvedemarken\';

n=0;
n=n+1;file_in{n}='HXB_E6789_R3_M15_IPproc7g_AutoProc_noiseGF_reproc030720.tx2';
n=n+1;file_in{n}='HXB_E6789_R4_IPproc7g_Res15_c_AutoProc.tx2';
n=n+1;file_in{n}='HXB_E6789_R5_7g_fObStd_M_15_AutoProc_noiseGF.tx2';

Mapp_cell=cell(n,1);
depths_cell=cell(n,1);
IP_flag_cell=cell(n,1);
t_center=cell(n,1);
current_cell=cell(n,1);
Rho_save=cell(n,1);
Res_save=cell(n,1);

% name_final='Hvede_manual_auto_deg4_normr1_v2.mat';
% load(name_final,'comp_dp','count_both_kept','count_both_removed','count_discr','auto_kept_not_manual','count_auto_kept_not_manual','manual_kept_notauto','comp_gates','count_gates_both_kept','count_gates_both_removed','count_gates_discr')


%% test the code
keep_all=1;
num_dataset=numel(file_in);
num_depth_ref = 2232;
my_voltage=cell(num_dataset,1);
time_curvature_max_all=cell(num_dataset);
fitresult_se=cell(num_dataset,num_depth_ref);
gof_se=cell(num_dataset,num_depth_ref);
output_se=cell(num_dataset,num_depth_ref);
se_fit_mapp=cell(num_dataset,num_depth_ref);
fitresult_se_2=cell(num_dataset,num_depth_ref);
gof_se_2=cell(num_dataset,num_depth_ref);
output_se_2=cell(num_dataset,num_depth_ref);
se_fit_mapp_2=cell(num_dataset,num_depth_ref);
coeff_alpha=NaN(num_dataset,num_depth_ref);
coeff_tauKWW=NaN(num_dataset,num_depth_ref);
coeff_beta=NaN(num_dataset,num_depth_ref);
id_not_flagged_EM=cell(num_dataset,num_depth_ref);
xdata0_excl_EM=cell(num_dataset,num_depth_ref);
ydata0_excl_EM=cell(num_dataset,num_depth_ref);
remaining_gates=cell(num_dataset,num_depth_ref);
remaining_m=cell(num_dataset,num_depth_ref);
val_clog_peak_all=NaN(num_dataset,num_depth_ref);
all_locs_all=NaN(num_dataset,num_depth_ref);
            
res_thres=0.1;
% rel_res_thres=0.02; %2%
rel_res_poly_thres=0.2;
lim_normr=1;
% lim_normr_end=0.5;
deg_pol=5;
min_peak_early=0.05;
min_peak_late=1;
max_time_peak_early=200;
plot_interm=0;
col=jet(num_dataset);

for k=2%1:num_dataset
    data=importdata(file_in{k});
    data=data.data;
    Ndata=size(data,1);
    Rho=data(:,22);
    Res=data(:,21);
    num_dp=Ndata;
    NGates=data(1,25);
    flagIP=data(:,end-8-NGates:end-9);
    delay=data(1,26+NGates);
    GateWidth=NaN(n,NGates);
    GateWidth(k,:) = data(1,26+NGates+1:26+2*NGates);
    GateStart=NaN(1,NGates);
    GateTime=NaN(1,NGates);
    GateStart(1)=delay;
    for i=2:NGates
        GateStart(i)=GateStart(i-1)+GateWidth(k,i-1);
    end
    for i=1:NGates
        GateTime(i)=GateStart(i)+0.5*GateWidth(k,i);
    end
    Current=data(:,end-7);    
    IP_Mapp = data(:,26:25+NGates);
    Ax=data(:,1);
    Bx=data(:,2);
    Mx=data(:,3);
    Nx=data(:,4);
    Az=data(:,5);
    Bz=data(:,6);
    Mz=data(:,8);
    Nz=data(:,4);

    current_cell{k}=Current;
    Mapp_cell{k}=IP_Mapp;
    IP_flag_cell{k}=flagIP;
    t_center{k}=GateTime;
    Rho_save{k}=Rho;
    Res_save{k}=Res;
    depths_cell{k}=[Ax Bx Mx Nx Az Bz Mz Nz];

    % calculate voltage
    my_voltage{k}=current_cell{k}.*Res_save{k};
    time_curvature_max=NaN(num_dp,1);
    time_cut_end=NaN(num_dp,1);
    time_curvature_max_save=zeros(num_dp,1);

    for i=304%1:num_dp  %ISL10: 1422%1418%1571%1495%1397%1433%1445%1502  %ISL12:  7%1383% %ISL11:728%1402% %ISL1:7%176%173%4%12%      
        if keep_all==1%manual_kept_notauto(k,i) %auto_kept_not_manual(k,i)%
%             f1 = figure(1);
%             clf(f1)
            clf;
            t = t_center{k};             % time
            m = Mapp_cell{k}(i,:);  % Chargeability
            negm = find(m<0);                 % remove negative m
            m(negm) = [];
            t(negm) = [];
            if numel(t)>5
                my_x=log(t);
                my_y=log(m);
                logt = linspace(min(my_x),max(my_x),100); % more data points on the fitted curve

                [p,S,mu] = polyfit(my_x,my_y,deg_pol);     % fitting the decay in log-log space
                logm_data = polyval(p,my_x,S,mu);  

                if plot_interm
                    my_time=exp(logt);
                    logm = polyval(p,logt,S,mu);           % evaluate logm on a dense logt
%                     hAx1 = axes;
%                     gca=hAx1;
                    yyaxis left

                    loglog(t_center{k}, abs(Mapp_cell{k}(i,:)),'or','MarkerFaceColor','r','MarkerSize',3);
                    hold on
                    loglog(t_center{k}, Mapp_cell{k}(i,:),'ob','MarkerFaceColor','b','MarkerSize',3)

                    xlim([1e0,1e4])
                    ylim([1e-1,1e4])
                    plot(my_time,exp(logm),'--b','LineWidth',0.5)
                    my_text=['S.normr=' num2str(S.normr)];
                    text(1e3,1e3,my_text);
                    hXLabel1=xlabel('Time (msec)');
                    hYLabel1=ylabel('Polarizability \eta (mV/V)');
                    hTitle1=title([file_in{k}(1:12) ' - DP' num2str(i)],'Interpreter','none');
%                     keyboard
                end

                while S.normr>lim_normr
                    rel_res_poly=(logm_data-my_y)./my_y;
                    excl_tmp = find(abs(rel_res_poly)>rel_res_poly_thres);
                    if isempty(excl_tmp)
                        break
                    else
                        newt=find(abs(rel_res_poly)<rel_res_poly_thres);
                        test_gate_early=find(newt<min(excl_tmp));
                        if ~isempty(test_gate_early)
                            if numel(test_gate_early)<4
                                newt(test_gate_early)=[]; % évite de garder des ealy gate si d'autres après ont été exclues
                            end
    %                             keyboard
                        end

                        m_poly=m(newt);
                        t_poly=t(newt);
                        if numel(t_poly)<=5
                            break
                        else
                            my_x_poly=log(t_poly);
                            my_y_poly=log(m_poly);
                            t=t_poly; % à enlever si besoin
                            m=m_poly;
                            my_x=my_x_poly;
                            my_y=my_y_poly;

                            logt = linspace(min(my_x_poly),max(my_x_poly),100); % more data points on the fitted curve
                            [p,S,mu] = polyfit(my_x_poly,my_y_poly,deg_pol);     % fitting the decay in log-log space
                            logm_data = polyval(p,my_x_poly,S,mu);  
                        end
                    end
                end

                rel_res_poly=(logm_data-my_y)./my_y;
                excl_tmp = find(abs(rel_res_poly)>rel_res_poly_thres);
                if isempty(excl_tmp)
                    m_poly=m;
                    t_poly=t;
                else
                    newt=find(abs(rel_res_poly)<rel_res_poly_thres);
                    test_gate_early=find(newt<min(excl_tmp));
                    if ~isempty(test_gate_early)
                        if numel(test_gate_early)<4
                            newt(test_gate_early)=[]; % évite de garder des ealy gate si d'autres après ont été exclues
                        end
                    end
                    m_poly=m(newt);
                    t_poly=t(newt);
                end        

                if numel(t_poly)>5
                    my_x_poly=log(t_poly);
                    my_y_poly=log(m_poly);
                    t=t_poly; % à enlever si besoin
                    m=m_poly;
                    my_x=my_x_poly;
                    my_y=my_y_poly;
                    logt = linspace(min(my_x_poly),max(my_x_poly),100); % more data points on the fitted curve
                    [p,S,mu] = polyfit(my_x_poly,my_y_poly,deg_pol);     % fitting the decay in log-log space
                    logm_data = polyval(p,my_x_poly,S,mu);  
                    logm = polyval(p,logt,S,mu);           % evaluate logm on a dense logt
                    clog = curL(logt,logm);
                    my_time=exp(logt);

                    TFmax = islocalmax(clog);
                    TFmin = islocalmin(clog);
                    locs_max=my_time(TFmax);
                    locs_min=my_time(TFmin);
                    val_clog_peak=[clog(TFmin),clog(TFmax)];
                    all_locs=[locs_min,locs_max];
                    locs=min(all_locs);
                    pos_locs=find(all_locs==locs);
                    locs_end=max(all_locs);
                    pos_locs_end=find(all_locs==locs_end);
                end


                if plot_interm
                    clf;
%                     hAx1 = axes;
%                     gca=hAx1;
                    yyaxis left

                    loglog(t_center{k}, abs(Mapp_cell{k}(i,:)),'or','MarkerFaceColor','r','MarkerSize',3);
                    hold on
                    loglog(t_center{k}, Mapp_cell{k}(i,:),'ob','MarkerFaceColor','b','MarkerSize',3)

                    xlim([1e0,1e4])
                    ylim([1e-1,1e4])
                    plot(my_time,exp(logm),'--b','LineWidth',0.5)
                    hXLabel1=xlabel('Time (msec)');
                    hYLabel1=ylabel('Polarizability \eta (mV/V)');
                    hTitle1=title([file_in{k}(1:12) ' - DP' num2str(i)],'Interpreter','none');
                    yyaxis right
                    plot(my_time,clog)
                    hYLabel2=ylabel('Curvature');
                    hold on

                    grid on
                    text(all_locs+.02,val_clog_peak,num2str((1:numel(val_clog_peak))'))
                    plot(my_time,clog,all_locs,val_clog_peak,'r*')

                    set(gca,'FontName','Times New Roman','FontSize',14);
                    set([hXLabel1,hYLabel1,hYLabel2,hTitle1],'FontName','Times New Roman','FontSize',14);

                    yyaxis right
                    hold off
                    yyaxis left
                    hold off
%                     keyboard
                end

                if isempty(locs)
                    warning(['no max found for file' file_in{k}(1:12) ' - DP ' num2str(i)])
                elseif numel(t)>5
                    val_clog_peak_all(k,i)=val_clog_peak(pos_locs);
                    all_locs_all(k,i)=locs;

                    if abs(val_clog_peak(pos_locs))<min_peak_early
                        time_curvature_max(i)=min(t);
                        ind_times_min = (t>time_curvature_max(i));
                    elseif locs>max_time_peak_early % to avoid removing something if no EM noise
                        time_curvature_max(i)=min(t);
                        ind_times_min = (t>time_curvature_max(i));                   
                    else
                        time_curvature_max(i)=locs;
                        ind_times_min = (t>time_curvature_max(i));
                        ind_tmp=find(ind_times_min);
                        if sum(ind_tmp)>0
                            ind_times_min(ind_tmp(1))=0; % to also cut the gate where the peak happens 
                        end
                    end   

                    if abs(val_clog_peak(pos_locs_end))<min_peak_late % check for weird shape at late times.
                        time_cut_end(i)=max(t);
                        ind_times_max = (t<time_cut_end(i));
                    elseif locs_end==locs % un seul pic, déjà compté comme early
                        time_cut_end(i)=max(t);
                        ind_times_max = (t<time_cut_end(i));
                    else
                        time_cut_end(i)=locs_end;
                        ind_times_max = (t<time_cut_end(i));
                        ind_tmp=find(ind_times_max);
                        ind_times_max(ind_tmp(end))=0; % to also cut the gate where the peak happens   
                    end            

                    ind_times_remaining = ind_times_min & ind_times_max;

                    remaining_gates{k,i}=t(ind_times_remaining); 
                    remaining_m{k,i}=m(ind_times_remaining);
                    remaining_gates_log=my_x(ind_times_remaining); % my_x is the log of time
                    remaining_m_log=my_y(ind_times_remaining);
                    if numel(remaining_gates{k,i})>6
                        [xData_exp, yData_exp] = prepareCurveData( remaining_gates_log, remaining_m_log);
                        ft = fittype( 'a1 - exp(c*x-b1)', 'independent', 'x', 'dependent', 'y' );
                        opts = fitoptions( 'Method', 'NonlinearLeastSquares');
                        opts.Display = 'Off';
                        opts.Robust = 'Off';
                        opts.Lower = [0 -Inf 0.1];
                        opts.StartPoint = [1 1 0.5];%[0.6 0.2 0.7];%[5 5 0.5];%
                        opts.Upper = [10 Inf 1];
                        % keyboard
                        [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                        coeff_alpha(k,i)=exp(fitresult_se{k,i}.a1);
                        coeff_tauKWW(k,i)=exp(fitresult_se{k,i}.b1/fitresult_se{k,i}.c);
                        coeff_beta(k,i)=fitresult_se{k,i}.c;
                        se_fit_mapp{k,i}=coeff_alpha(k,i)*exp(-(remaining_gates{k,i}/coeff_tauKWW(k,i)).^coeff_beta(k,i));

                        if gof_se{k, i}.rsquare<0.9
                            warning(['Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:12) '; try new starting point'])
                            opts.StartPoint = [5 5 0.5];
                            [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                            coeff_alpha(k,i)=exp(fitresult_se{k,i}.a1);
                            coeff_tauKWW(k,i)=exp(fitresult_se{k,i}.b1/fitresult_se{k,i}.c);
                            coeff_beta(k,i)=fitresult_se{k,i}.c;
                            se_fit_mapp{k,i}=coeff_alpha(k,i)*exp(-(remaining_gates{k,i}/coeff_tauKWW(k,i)).^coeff_beta(k,i));
                            if gof_se{k, i}.rsquare<0.9
                                warning(['(again) Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:12) '; try new starting point'])
                                opts.StartPoint = [0.6 0.2 0.7];
                                [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                                coeff_alpha(k,i)=exp(fitresult_se{k,i}.a1);
                                coeff_tauKWW(k,i)=exp(fitresult_se{k,i}.b1/fitresult_se{k,i}.c);
                                coeff_beta(k,i)=fitresult_se{k,i}.c;
                                se_fit_mapp{k,i}=coeff_alpha(k,i)*exp(-(remaining_gates{k,i}/coeff_tauKWW(k,i)).^coeff_beta(k,i));
                            end
                        end

                        if gof_se{k, i}.rsquare>0.9
                            time_curvature_max_save(i)=1;
                        end

                        if gof_se{k, i}.rsquare>0.99                        
                            time_curvature_max_save(i)=1; % ne sert à rien mais OK
                        elseif isempty(fitresult_se{k,i})
                            warning(['no SE fit for DP =' num2str(i) 'm for file' file_in{k}  '- with EM'])
                            id_not_flagged_EM{k,i}=zeros(1,NGates); % 1 if not flagged
                            EM_ok=0;
                        else
                            % EM_ok=1;
                            residuals_Krafla_EM = output_se{k, i}.residuals; 
    %                         relative_residuals_Krafla_EM = residuals_Krafla_EM./se_fit_mapp{k,i}';
    %                         points_to_keep_EM=find(abs(relative_residuals_Krafla_EM)<rel_res_thres);

                            points_to_keep_EM=find(abs(residuals_Krafla_EM)<res_thres);

                            xdata0_excl_EM{k,i}=remaining_gates{k,i}(points_to_keep_EM);
                            ydata0_excl_EM{k,i}=remaining_m{k,i}(points_to_keep_EM);
                            id_not_flagged_EM{k,i}=ismember(t_center{k}(1:end),xdata0_excl_EM{k, i}); % 1 if not flagged                       

                            if numel(xdata0_excl_EM{k,i})>6
                                xdata0_excl_EM_log=log(xdata0_excl_EM{k,i});
                                ydata0_excl_EM_log=log(ydata0_excl_EM{k,i});
                                [xData_exp, yData_exp] = prepareCurveData( xdata0_excl_EM_log, ydata0_excl_EM_log);
                                ft = fittype( 'a1 - exp(c*x-b1)', 'independent', 'x', 'dependent', 'y' );
                                opts = fitoptions( 'Method', 'NonlinearLeastSquares');
                                opts.Display = 'Off';
                                opts.Robust = 'Off';
                                opts.Lower = [0 -Inf 0.1];
                                opts.StartPoint = [1 1 0.5];%[0.6 0.2 0.7];%[5 5 0.5];%
                                opts.Upper = [10 Inf 1];
                                % keyboard
                                [fitresult_se_2{k,i}, gof_se_2{k,i}, output_se_2{k,i}] = fit( xData_exp, yData_exp, ft, opts );

                                if gof_se_2{k, i}.rsquare<0.9
                                    warning(['Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:12) '; try new starting point'])
                                    opts.StartPoint = [5 5 0.5];
                                    [fitresult_se_2{k,i}, gof_se_2{k,i}, output_se_2{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                                    if gof_se_2{k, i}.rsquare<0.9
                                        warning(['(again) Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:12) '; try new starting point'])
                                        opts.StartPoint = [0.8 0.6 0.3];
                                        [fitresult_se_2{k,i}, gof_se_2{k,i}, output_se_2{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                                    end
                                end
                                coeff_alpha(k,i)=exp(fitresult_se_2{k,i}.a1);
                                coeff_tauKWW(k,i)=exp(fitresult_se_2{k,i}.b1/fitresult_se_2{k,i}.c);
                                coeff_beta(k,i)=fitresult_se_2{k,i}.c;       
                                se_fit_mapp_2{k,i}=coeff_alpha(k,i)*exp(-(xdata0_excl_EM{k,i}/coeff_tauKWW(k,i)).^coeff_beta(k,i));
                                if gof_se_2{k, i}.rsquare>0.9
                                    time_curvature_max_save(i)=1;
                                end
                            end
                        end
                    end
                else
                    warning(['not enough data left for polyfit' file_in{k}(1:12) ' - DP ' num2str(i)])
                end

                if sum(Mapp_cell{k}(i,:))==0
                    warning(['no IP data for' file_in{k}(1:12) ' - DP ' num2str(i)])
                else
                    clf;

    %             if ismember(i,[1675,1417,1422,1413,1673,1408,1383,1669,1656,1410,1402])
%                     f1 = figure(1);
%                     hAx1 = axes;
%                     gca=hAx1;
%                     gcf=f1;
                    yyaxis left
                    % loglog(hAx1,t,m,'ok','MarkerSize',8)

                    loglog(t_center{k}, abs(Mapp_cell{k}(i,:)),'or','MarkerFaceColor','r','MarkerSize',3)
                    hold on
                    loglog(t_center{k}, Mapp_cell{k}(i,:),'ob','MarkerFaceColor','b','MarkerSize',3)

                    xlim([1e0,1e4])
                    ylim([1e-1,1e4])
                    plot(my_time,exp(logm),'--b','LineWidth',0.5)
                    hXLabel1=xlabel('Time (msec)');
                    hYLabel1=ylabel('Polarizability \eta (mV/V)');

                    plot(remaining_gates{k,i},remaining_m{k,i},'o','Color','g','MarkerFaceColor',[0 1 0],'MarkerSize',6);

                    if ~isempty(gof_se{k, i})
                        if numel(xdata0_excl_EM{k,i})>6
                            plot(remaining_gates{k,i},remaining_m{k,i},'x','Color','k','MarkerSize',8);
                            plot(xdata0_excl_EM{k,i},ydata0_excl_EM{k,i},'o','Color','k','MarkerFaceColor','k','MarkerSize',6);
                            plot(xdata0_excl_EM{k,i},se_fit_mapp_2{k,i},'-','Color','k','LineWidth',2);
                            my_txt=sprintf('Cut at %.2f ms \n R^2_{SE2} = %.4f \n R^2_{SE1} = %.4f',time_curvature_max(i), gof_se_2{k, i}.rsquare, gof_se{k, i}.rsquare);
                            hText=text(500,2000,my_txt);
                            set(hText,'FontName','Times New Roman','FontSize',12);
                        elseif gof_se{k, i}.rsquare>0.9                       
                            plot(remaining_gates{k,i},se_fit_mapp{k,i},'-','Color','k','LineWidth',2)
                            my_txt=sprintf('Cut at %.2f ms \n R^2_{SE1} = %.4f',time_curvature_max(i), gof_se{k, i}.rsquare);
                            hText=text(500,2000,my_txt);
                            set(hText,'FontName','Times New Roman','FontSize',12);
                        end
                    end

                    hTitle1=title([file_in{k}(1:12) ' - DP' num2str(i)],'Interpreter','none');
                    yyaxis right
                    plot(my_time,clog)
                    hYLabel2=ylabel('Curvature');
                    hold on

                    grid on
                    text(all_locs+.02,val_clog_peak,num2str((1:numel(val_clog_peak))'))
                    plot(my_time,clog,all_locs,val_clog_peak,'r*')

                    set(gca,'FontName','Times New Roman','FontSize',14);
                    set([hXLabel1,hYLabel1,hYLabel2,hTitle1],'FontName','Times New Roman','FontSize',14);

                    yyaxis right
                    hold off
                    yyaxis left
                    hold off

    %                 name_out=['png\Krafla_SEfit_EM_deg4_normr1\' file_in{k}(1:12) '_DP' num2str(i) '_SEfit_EM_deg4_v3.png'];
%                     name_out=['png\Hvede_SEfit_EM_deg5_manual_kept_not_autodeg4\' file_in{k}(1:12) '_DP' num2str(i) '_SEfit_EM_deg5.png'];
%                     name_out=['png\Hvede_SEfit_EM_deg5\' file_in{k}(1:12) '_DP' num2str(i) '_SEfit_EM_deg5.png'];
                    name_out=['pdf\' file_in{k}(1:12) '_DP' num2str(i) '_SEfit_EM_deg5.pdf'];

                    exportgraphics(gcf,name_out,'Resolution',300)
    %                 keyboard
                end
            else
                warning(['Not enough gates for file' file_in{k}(1:12) ' - DP ' num2str(i)])
            end
%             keyboard
        end
    end
    time_curvature_max_all{k}=[time_curvature_max,time_curvature_max_save];
end


%% plot time vs voltage
f2 = figure(2);
hAx2 = axes;
for k=1:num_dataset
    if ~isempty(my_voltage{k})
        plot(hAx2,my_voltage{k}(time_curvature_max_all{k}(:,2)==1),time_curvature_max_all{k}(time_curvature_max_all{k}(:,2)==1,1),'ok','MarkerFaceColor',col(k,:))
        hold(hAx2,'on')
    end
end
set(hAx2, 'XScale', 'log')
set(hAx2, 'YScale', 'log')
% xlim([1e-3,1e2])
xlim([1e-4,1e1])
ylim([1e0,1e2])
grid on;
hXLabel=xlabel('Voltage (V)');
% hXLabel=xlabel('Resistance (\Omega)');
hYLabel=ylabel('Time of curvature max (ms)');
hLegend = legend(file_in, 'Interpreter','none'); 
gca=hAx2;
set(gca,'FontName','Times New Roman','FontSize',14);
set([hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
gcf=f2;
name_out2='png\Hvede_EM_voltage_v9.png';
exportgraphics(gcf,name_out2,'Resolution',300)
% clf(f2);


%% histograms for curl max
% for histogram of peaks val_clog_peak_all(k,i)
% for histograms of times all_locs_all(k,i)

%% save mat file
name_final='Hvede_SEfit_log_curvatureMax_Min_End_deg5.mat';
save(name_final,'num_depth_ref','remaining_gates','remaining_m','xdata0_excl_EM','ydata0_excl_EM','val_clog_peak_all','all_locs_all','time_curvature_max_all','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2','coeff_alpha','coeff_tauKWW', 'coeff_beta','Mapp_cell','IP_flag_cell','t_center','current_cell','my_voltage','Rho_save','Res_save','file_in')

