close all
clear all
clc



load  ukb\wmh\ukb_subsets\Cognitive_function_data.mat

Cognitive_function.Cognitive_names={'Pairs matching'
    'Numeric memory'
    'Trail making'
    'Matrix pattern completion'
    'Fluid intelligence'
    'Prospective memory'
    'Reaction time'
    'Tower rearranging'
    'Symbol digit substitution'};

Cognitive_names=Cognitive_function.Cognitive_names;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wmh_cov_info_0=readtable('ukb\wmh/wmh_cov_info.csv');%,'WriteRowNames',true

ind2=sum(isnan(wmh_cov_info_0(:,[2:end]).Variables'))==0;
wmh_cov_info_0=wmh_cov_info_0(ind2,:);


load  ukb\wmh/wmh_sur_cov.mat

heart_disease=wmh_sur_cov(:,[end-3:end]).Variables;
HD_ind=sum(heart_disease')>0;

wmh_sur_cov=wmh_sur_cov(:,[1:end-3]);
wmh_sur_cov.heart_disease=HD_ind';

load  ukb\wmh/wmh_cov_data.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pred_wmh=readtable('s5_Pred.csv');


[final_eid,~]= intersect(intersect(intersect(intersect(wmh_sur_cov.eid,Cognitive_function_Instance2_data.eid),pred_wmh.eid),cov_info_2.eid),wmh_cov_info_0.eid);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[~,ind2]= intersect(Cognitive_function_Instance2_data.eid,final_eid);
Cognitive_function_Instance2_data1=Cognitive_function_Instance2_data(ind2,:);


 [~,ind2]= intersect(wmh_sur_cov.eid,final_eid);
 wmh_sur_cov=wmh_sur_cov(ind2,:);


[~,ind2]= intersect(pred_wmh.eid,final_eid);
pred_wmh=pred_wmh(ind2,:);


[~,ind2]= intersect(wmh_cov_info_0.eid,final_eid);
wmh_cov_info_0=wmh_cov_info_0(ind2,:);

[~,ind2]= intersect(cov_info_2.eid,final_eid);
cov_info_2=cov_info_2(ind2,:);




data_eid=[wmh_sur_cov.eid,Cognitive_function_Instance2_data1.eid,wmh_cov_info_0.eid,pred_wmh.eid];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


site_num = unique(wmh_cov_info_0.centre);
num_each_size=[];
for i=1:length(site_num)
    num_each_size(i) =  sum(wmh_cov_info_0.centre==site_num(i));
end
bl_site = zeros(length(wmh_cov_info_0.centre),sum(num_each_size>100));
for i=1:length(site_num)
    if num_each_size(i)>100
        bl_site(wmh_cov_info_0.centre==site_num(i),i) = 1;
    end
end
index = sum(bl_site);
bl_site(:,index==0) = [];
bl_site(:,end) = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wmh_cov_info_0.Ethnic2=wmh_cov_info_0.Ethnic==1001;
bl_Ethnic=wmh_cov_info_0.Ethnic==1001;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bl_Medical_covinfo=wmh_sur_cov(:,[2:end]).Variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
instance2_site_num = unique(cov_info_2.x54_2_0);
num_site2_size=[];
for i=1:length(instance2_site_num)
    num_site2_size(i) =  sum(cov_info_2.x54_2_0==instance2_site_num(i));
end
instance2_site = zeros(length(cov_info_2.x54_2_0),sum(num_site2_size>100));
for i=1:length(instance2_site_num)
    if num_site2_size(i)>100
        instance2_site(cov_info_2.x54_2_0==instance2_site_num(i),i) = 1;
    end
end
index = sum(instance2_site);
instance2_site(:,index==0) = [];
instance2_site(:,end) = [];



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
instance2_duration=nan(size(cov_info_2,1),1);
for i =1: size(cov_info_2,1)

    bl_date=cov_info_bl(i,:).x53_0_0;
    instance2_date=cov_info_2(i,:).x53_2_0;

    t_between= between(bl_date, instance2_date,'months');
    instance2_duration(i,1)=split(t_between,'mo');
end

instance2_age=wmh_cov_info_0.age+instance2_duration./12;


%%

cmap=colormap(lines);

instance2_cov_info=[wmh_cov_info_0.sex,zscore(instance2_age), bl_Ethnic,bl_site,instance2_site,bl_Medical_covinfo];%,



predict_wmh_corr_Cognitive_results=nan(size(Cognitive_names,1),2);
corr_rval=nan(size(Cognitive_names,1),1);
corr_pval=nan(size(Cognitive_names,1),1);

ValueType=[];
for i=1:size(Cognitive_names,1)

    tic
 
    y1=Cognitive_function_Instance2_data1(:,i+1).Variables;
    ind2=~isnan(y1);

    y=y1(ind2,:);
    x=pred_wmh(ind2,:).y_pred;

    cov_info2=instance2_cov_info(ind2,:);
    index = sum(cov_info2);
    cov_info2(:,index==0) = [];
    z=cov_info2;

    if i==3
           % TF = isoutlier(y,"quartiles");

            TF =y>1000;
        
            x=x(TF==0,:);
            y=y(TF==0,:);
            z=z(TF==0,:);
    end
    nn(i)=length(x);

    %x=zscore(x)
    %y=zscore(y)
    [corr_rval(i,1),corr_pval(i,1)]= partialcorr(x,y,z,'type' ,'Spearman');% );%

    corr_pval(i,1)=corr_pval(i,1) *9;


    figure(i)
    %subplot(2,3,k)
    sz=1;
    scatter(x,y,sz,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',1)  %'filled'
    h1 = lsline;
    h1.Color = 'r';
    h1.LineWidth = 2.5


    aa= (max(y)-min(y))/10;

    ylim([min(y)-aa max(y)+aa*1]);
    ylim([min(y)-aa max(y)+aa*1.2]);

    xlim([min(x)-0.2 max(x)+0.2]);


  
    set(gca,'fontname','times')  % Set it to times
    set(gcf,'position',[300 300 450 400])
    set(gca,'TickDir' , 'out', 'TickLength', [0.02 0.25])
    set(gca,'linewidth',1)
    box off


    str=sprintf(Cognitive_names{i});
    ylabel(str,'Fontsize',14)
    xlabel(['Predicted WMH'],'Fontsize',14)


    xt=min(x)+0.5;
    yt=max(y)+ aa/1.5;

    if corr_pval(i)> 0.001
        text(xt,yt,{['\rmR = ' sprintf('%4.2f', corr_rval(i)),...
            '   \itP\rm = ', sprintf('%0.2g',corr_pval(i))]},'Fontsize',14,'fontname','times')
    else
        text(xt,yt,{['\rmR = ' sprintf('%4.2f', corr_rval(i)), ...
            '   \itP\rm < 0.001 ']},'Fontsize',14,'fontname','times')
    end
%  ' \rmN = ', num2str(nn(i))]
    print(figure(i), '-dtiff','-r600', ['figs/predicted_wmh_corr_',Cognitive_names{i},'.tif'])

    toc

end


predict_wmh_corr_Cognitive_results=array2table([corr_rval,corr_pval],"VariableNames",{'rvalue','raw_p'},'RowNames', Cognitive_names);
predict_wmh_corr_Cognitive_results.Cognitive_names=Cognitive_names;

