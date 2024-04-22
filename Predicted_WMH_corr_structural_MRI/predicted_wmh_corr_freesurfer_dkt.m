
close all
clear all
clc



cd  ukb\wmh\wmh_corr_structural_mri

load ukb\t1_mri\ukb_data\Freesurfer_DKT_data.mat
Freesurfer_DKT_names=Freesurfer_DKT.Description;


x=Freesurfer_DKT_bl_data.x26528_2_0;
y=Freesurfer_DKT_bl_data.x26562_2_0;

ind2=sum(~isnan([x,y]'))>1;
[r1,p1]=corr(x(ind2),y(ind2))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %
ukb_wmh=readtable('ukb\wmh\YJ/s5_Pred.csv');
%ind2=sum(~isnan(ukb_wmh(:,[2:end]).Variables'))>0;
ind2=~isnan(ukb_wmh.y_pred)>0;
ukb_wmh=ukb_wmh(ind2,:);


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

[final_eid,~]= intersect(intersect(intersect(intersect(ukb_wmh.eid,Freesurfer_DKT_bl_data.eid),wmh_cov_info_0.eid),wmh_sur_cov.eid),cov_info_2.eid);

%%
[~,ind2]= intersect(Freesurfer_DKT_bl_data.eid,final_eid);
Freesurfer_DKT_bl_data1=Freesurfer_DKT_bl_data(ind2,:);


[~,ind2]= intersect(ukb_wmh.eid,final_eid);
ukb_wmh=ukb_wmh(ind2,:);



[~,ind2]= intersect(wmh_cov_info_0.eid,final_eid);
wmh_cov_info_0=wmh_cov_info_0(ind2,:);

[~,ind2]= intersect(wmh_sur_cov.eid,final_eid);
wmh_sur_cov=wmh_sur_cov(ind2,:);

[~,ind2]= intersect(cov_info_2.eid,final_eid);
cov_info_2=cov_info_2(ind2,:);


data_eid=[ukb_wmh.eid,Freesurfer_DKT_bl_data1.eid,wmh_cov_info_0.eid];

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



instance2_duration=nan(size(cov_info_2,1),1);
for i =1: size(cov_info_2,1)

    bl_date=cov_info_bl(i,:).x53_0_0;
    instance2_date=cov_info_2(i,:).x53_2_0;

    t_between= between(bl_date, instance2_date,'months');
    instance2_duration(i,1)=split(t_between,'mo');
end

%instance2_duration=zscore(instance2_duration./12);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%   

instance2_age=wmh_cov_info_0.age+instance2_duration./12;


% 26518	'Volume of TotalGray (whole brain)'
instance2_cov_info=[wmh_cov_info_0.sex,zscore(instance2_age), bl_Ethnic,bl_site,instance2_site,bl_Medical_covinfo,Freesurfer_DKT_bl_data1.x26521_2_0];%,


wmh_corr_dkt_results=nan(size(Freesurfer_DKT_names,1),2);
corr_rval=nan(size(Freesurfer_DKT_names,1),1);
corr_pval=nan(size(Freesurfer_DKT_names,1),1);

ValueType=[];
for i=1:size(Freesurfer_DKT_names,1)

    tic

    x1=Freesurfer_DKT_bl_data1(:,i+1).Variables;
    ind2=~isnan(x1);

    x=x1(ind2,:);
    y=ukb_wmh(ind2,:).y_pred;
    z=instance2_cov_info(ind2,:);

    [corr_rval(i,1),corr_pval(i,1)]= partialcorr(x,y,z);% );%

    toc

end



wmh_corr_dkt_results=array2table([corr_rval,corr_pval],"VariableNames",{'rvalue','raw_p'},'RowNames', Freesurfer_DKT_names);
wmh_corr_dkt_results.Freesurfer_DKT_names=Freesurfer_DKT_names;




writetable(wmh_corr_dkt_results,'ukb/wmh/wmh_corr_structural_mri/results/predict_wmh_corr_dkt_results.csv');








