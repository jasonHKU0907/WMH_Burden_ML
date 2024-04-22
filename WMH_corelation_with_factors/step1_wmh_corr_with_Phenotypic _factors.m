
close all
clear all
clc

addpath('ukb\wmh')

cd  ukb\wmh\Cross_section_correlation
load ukb\wmh/WMH_Phenotypic_bl_data.mat


%Phenotypic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Feature_Dict=readtable('ukb\UKB_YJ\Feature_Dict_2.csv');

Feature_Dict2=[];
for i=1:length(Phenotypic_bl_names.FieldID)
    [~, ind2] = intersect(Feature_Dict.FieldID,Phenotypic_bl_names(i,:).FieldID);
    Feature_Dict2=[Feature_Dict2;Feature_Dict(ind2,:)];
end

Phenotypic_bl_names.FieldID3=Feature_Dict2.FieldID;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ukb_wmh=readtable('../WMH_ukb_predict.csv');

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


[final_eid,~]= intersect(intersect(intersect(ukb_wmh.eid,Phenotypic_bl_data.eid),wmh_cov_info_0.eid),wmh_sur_cov.eid);

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,ind2]= intersect(Phenotypic_bl_data.eid,final_eid);
Phenotypic_bl_data1=Phenotypic_bl_data(ind2,:);


[~,ind2]= intersect(ukb_wmh.eid,final_eid);
ukb_wmh=ukb_wmh(ind2,:);


[~,ind2]= intersect(wmh_cov_info_0.eid,final_eid);
wmh_cov_info_0=wmh_cov_info_0(ind2,:);


[~,ind2]= intersect(wmh_sur_cov.eid,final_eid);
wmh_sur_cov=wmh_sur_cov(ind2,:);


[~,ind2]= intersect(cov_info_2.eid,final_eid);
cov_info_2=cov_info_2(ind2,:);



data_eid=[ukb_wmh.eid,Phenotypic_bl_data1.eid,wmh_cov_info_0.eid];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

instance2_duration=instance2_duration./12;

a_m1=round(median(instance2_duration),1);
a_Q1 = round(quantile(instance2_duration, 0.25),1); 
a_Q3 = round(quantile(instance2_duration, 0.75),1); 

%instance2_duration=zscore(instance2_duration./12);


sum(~isnan(Phenotypic_bl_data1.Variables))/38219

%%   % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wmh_log2=real(log2(ukb_wmh.x25781_2_0));


instance2_cov_info=[wmh_cov_info_0.sex,zscore(wmh_cov_info_0.age), bl_Ethnic,instance2_site,instance2_duration,wmh_sur_cov(:,[2:end]).Variables];%,

index = sum(instance2_cov_info);
instance2_cov_info(:,index==0) = [];

[T,Res,total_wmh_res,beta]=BWAS_Tregression(instance2_cov_info,wmh_log2);



bl_cov_info=[wmh_cov_info_0.sex,zscore(wmh_cov_info_0.age), bl_Ethnic,bl_site,wmh_sur_cov(:,[2:end]).Variables];%,bl_Medical_covinfo



wmh_corr_Phenotypic_results=nan(size(Phenotypic_bl_names,1),2);
corr_rval=nan(size(Phenotypic_bl_names,1),1);
corr_pval=nan(size(Phenotypic_bl_names,1),1);

ValueType=[];
eid_num =nan(size(Phenotypic_bl_names,1),1);
trait_mean=[];
trait_sd=[];

for i=1:size(Phenotypic_bl_names,1)

    tic
    ValueType{i,1}=Feature_Dict2(i,:).ValueType;

    x1=Phenotypic_bl_data1(:,i+1).Variables;
    ind2=~isnan(x1);



    bl_cov_info2=bl_cov_info(ind2,:);
    index = sum(bl_cov_info2);
    bl_cov_info2(:,index==0) = [];

    if strcmp(ValueType{i,1},'Continuous')| strcmp(ValueType{i,1},'Integer')

        [T,Res,x1_res,beta]=BWAS_Tregression(bl_cov_info2,x1(ind2));
        x=x1_res;
    else
        x=x1(ind2,:);
    end

    y=total_wmh_res(ind2,:);

    [corr_rval(i,1),corr_pval(i,1)]= corr(x,y,'type' ,'Spearman');% );%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eid_num(i,1)=size(y,1);

    x1=Phenotypic_bl_data1(:,i+1).Variables;
    ind2=~isnan(x1);
    if strcmp(ValueType{i,1},'Continuous')  | strcmp(ValueType{i,1},'Integer')

        t_mean=round(nanmean(x1),3);
        trait_mean{i,1}=num2str(t_mean);

        t_std=round(nanstd(x1),3);
        trait_sd{i,1}=num2str(t_std);

    elseif max(x1)==1;

        t_mean=sum(x1==1);
        trait_mean{i,1}=num2str(t_mean);

        t_std=round(sum(x1==1)./eid_num(i,1),3) *100;
        trait_sd{i,1}=strcat([num2str(t_std),'%']);

    elseif max(x1)>1;
        t2_mean=[];
        t2_std=[];
        for n=1:max(x1)
            t_mean=sum(x1==n);
            t2_mean=strcat([t2_mean, num2str(t_mean),';']);

            t_std=round(sum(x1==n)./eid_num(i,1),3) *100;
            t2_std=strcat([t2_std, num2str(t_std),'%;']);
        end
        trait_mean{i,1}=t2_mean;
        trait_sd{i,1}=t2_std;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc

end



%p_fdr=mafdr(corr_pval,'BHFDR', true);

wmh_corr_Phenotypic_results=array2table([corr_rval,corr_pval],"VariableNames",{'rvalue','raw_p'},'RowNames', Phenotypic_bl_names.Description);
wmh_corr_Phenotypic_results.Phenotypic_bl_names=Phenotypic_bl_names.Description;
wmh_corr_Phenotypic_results.category=Phenotypic_bl_names.category;


wmh_corr_Phenotypic_table=wmh_corr_Phenotypic_results;

[~,ind2]=sort(wmh_corr_Phenotypic_results.category);
wmh_corr_Phenotypic_results=wmh_corr_Phenotypic_results(ind2,:);

writetable(wmh_corr_Phenotypic_results,'ukb/wmh/Cross_section_correlation/results/wmh_corr_Phenotypic_results.csv');


