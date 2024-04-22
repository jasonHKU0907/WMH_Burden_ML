
close all
clear all
clc

addpath('ukb\wmh')

cd  ukb\wmh\wmh_corr_structural_mri
load ukb\wmh/WMH_Phenotypic_bl_data.mat

Phenotypic_feildid=Phenotypic_bl_data(:,1);
Phenotypic_bl_data=Phenotypic_bl_data(:,[2:end]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Feature_Dict=readtable('ukb\UKB_YJ\Feature_Dict_2.csv');

Feature_Dict2=[];
for i=1:length(Phenotypic_bl_names.FieldID)
    [~, ind2] = intersect(Feature_Dict.FieldID,Phenotypic_bl_names(i,:).FieldID);
    Feature_Dict2=[Feature_Dict2;Feature_Dict(ind2,:)];
end

Phenotypic_bl_names.FieldID3=Feature_Dict2.FieldID;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wmh_corr_Phenotypic_results=readtable('D:\华为家庭存储\I-盘/ukb/wmh/Cross_section_correlation/results/wmh_corr_Phenotypic_results.csv');


[~,ind2]= intersect(wmh_corr_Phenotypic_results.Phenotypic_bl_names,Phenotypic_bl_names.Description);
wmh_corr_Phenotypic_results=wmh_corr_Phenotypic_results(ind2,:);


[~,ind2]= intersect(Phenotypic_bl_names.Description,wmh_corr_Phenotypic_results.Phenotypic_bl_names);
Phenotypic_bl_names=Phenotypic_bl_names(ind2,:);
Phenotypic_bl_data=Phenotypic_bl_data(:,ind2);


Phenotypic_bl_names.Phenotypic_bl_names3=wmh_corr_Phenotypic_results.Phenotypic_bl_names;
Phenotypic_bl_names.Phenotypic_bl_names4=Phenotypic_bl_data.Properties.VariableNames';

% %
ind2=wmh_corr_Phenotypic_results.raw_p<0.05/205;
Phenotypic_bl_data=Phenotypic_bl_data(:,ind2);
Phenotypic_bl_names=Phenotypic_bl_names(ind2,:);

Phenotypic_bl_data=[Phenotypic_feildid,Phenotypic_bl_data];



%%
category=unique(Phenotypic_bl_names.category);

num=[];
for i=1:length(category)

    ind2=strcmp(Phenotypic_bl_names.category,category{i});
    tt=Phenotypic_bl_names(ind2,:);

    num(i,1)=sum(ind2);

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load ukb\wmh\ukb_subsets/Diffusion_MRI_data.mat
Diffusion_MRI_names=Diffusion_MRI.Description;

dti_category=[];
side=[];
roi_names=[];

for k=1:length(Diffusion_MRI_names)
    d_name=split(Diffusion_MRI_names{k},' ');

    dti_category{k,1}=d_name{2};
    side{k,1}=d_name{end};

    d2_names=[];
    for n=1:length(d_name)-5;
        d2_names=strcat([d2_names,d_name{4+n},'_']);
    end

    d2_names=d2_names(1:end-1);
    if strcmp(d_name{end},'(right)')
        d2_names=strcat('R_',d2_names);
    elseif strcmp(d_name{end},'(left)')
        d2_names=strcat('L_',d2_names);

    elseif strcmp(d_name{end},'major')
        d2_names=strcat(d2_names,'_major');
    elseif strcmp(d_name{end},'minor')
        d2_names=strcat(d2_names,'_minor');
        %peduncle
    elseif strcmp(d_name{end},'peduncle')
        d2_names=strcat(d2_names,'_peduncle');
    end

    d3_names=strcat([upper(d2_names(1)),d2_names(2:end)]);
    roi_names{k,1}=d3_names;
end

Diffusion_MRI.dti_category=dti_category;
Diffusion_MRI.roi_names=roi_names;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

[final_eid,~]= intersect(intersect(intersect(intersect(Diffusion_MRI_bl_data.eid,Phenotypic_bl_data.eid),wmh_cov_info_0.eid),cov_info_2.eid),wmh_sur_cov.eid);

%%
[~,ind2]= intersect(Phenotypic_bl_data.eid,final_eid);
Phenotypic_bl_data1=Phenotypic_bl_data(ind2,:);


[~,ind2]= intersect(Diffusion_MRI_bl_data.eid,final_eid);
Diffusion_MRI_bl_data=Diffusion_MRI_bl_data(ind2,:);
%save('ukb_subsets/ukb673929','ukb673929', '-v7.3');


[~,ind2]= intersect(wmh_cov_info_0.eid,final_eid);
wmh_cov_info_0=wmh_cov_info_0(ind2,:);


[~,ind2]= intersect(wmh_sur_cov.eid,final_eid);
wmh_sur_cov=wmh_sur_cov(ind2,:);


[~,ind2]= intersect(cov_info_2.eid,final_eid);
cov_info_2=cov_info_2(ind2,:);



data_eid=[Diffusion_MRI_bl_data.eid,Phenotypic_bl_data1.eid,wmh_cov_info_0.eid];


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
instance2_duration=nan(size(cov_info_2,1),1);
for i =1: size(cov_info_2,1)

    bl_date=cov_info_bl(i,:).x53_0_0;
    instance2_date=cov_info_2(i,:).x53_2_0;

    t_between= between(bl_date, instance2_date,'months');
    instance2_duration(i,1)=split(t_between,'mo');
end

instance2_duration=zscore(instance2_duration./12);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sum(~isnan(Phenotypic_bl_data1.Variables))/38976;

%%




bl_cov_info=[wmh_cov_info_0.sex,zscore(wmh_cov_info_0.age), bl_Ethnic,bl_site,wmh_sur_cov(:,[2:end]).Variables];%,bl_Medical_covinfo
instance2_cov_info=[wmh_cov_info_0.sex,zscore(wmh_cov_info_0.age), bl_Ethnic,instance2_site,instance2_duration,wmh_sur_cov(:,[2:end]).Variables];%,


corr_rval=nan(size(Phenotypic_bl_names,1),size(Diffusion_MRI_bl_data,2)-1);
corr_pval=nan(size(Phenotypic_bl_names,1),size(Diffusion_MRI_bl_data,2)-1);


for k=1:size(Diffusion_MRI_bl_data,2)-1

    tic
    y1_dti=Diffusion_MRI_bl_data(:,k+1).Variables;


    ValueType=[];
    for i=1:size(Phenotypic_bl_names,1)


        ValueType{i,1}=Feature_Dict2(i,:).ValueType;
        x1=Phenotypic_bl_data1(:,i+1).Variables;

        ind2=~isnan([x1,y1_dti]);
        ind2=sum(ind2')>1;

        bl_cov_info2=bl_cov_info(ind2,:);
        index = sum(bl_cov_info2);
        bl_cov_info2(:,index==0) = [];

        if strcmp(ValueType{i,1},'Continuous') | strcmp(ValueType{i,1},'Integer')

            [T,Res,x1_res,beta]=BWAS_Tregression(bl_cov_info2,x1(ind2));
            x=x1_res;
        else
            x=x1(ind2,:);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %instance2_cov_info2=instance2_cov_info(ind2,:);
        [T,Res,y1_dti_res,beta]=BWAS_Tregression(instance2_cov_info(ind2,:),y1_dti(ind2,:));

        y=y1_dti_res;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [corr_rval(i,k),corr_pval(i,k)]= corr(x,y,'type' ,'Spearman');% );%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end



    toc
end
%p_fdr=mafdr(corr_pval,'BHFDR', true);



%%
Phenotypic_dti_corr_rval=array2table(corr_rval,"VariableNames",Diffusion_MRI_bl_data(:,[2:end]).Properties.VariableNames,'RowNames', Phenotypic_bl_names.Description);
Phenotypic_dti_corr_pval=array2table(corr_pval,"VariableNames",Diffusion_MRI_bl_data(:,[2:end]).Properties.VariableNames,'RowNames', Phenotypic_bl_names.Description);




%%% FA
ind2=strcmp(Diffusion_MRI.dti_category,'FA');
FA_names=Diffusion_MRI(ind2,:);

FA_rvalues=Phenotypic_dti_corr_rval(:,ind2);
FA_pvalues=Phenotypic_dti_corr_pval(:,ind2);

FA_names.FieldID2=FA_rvalues.Properties.VariableNames';

FA_rvalues.Properties.VariableNames=FA_names.roi_names;
FA_pvalues.Properties.VariableNames=FA_names.roi_names;

writetable(FA_rvalues,'ukb\wmh\wmh_corr_structural_mri/results/FA_rvalues.csv','WriteRowNames',true);
writetable(FA_pvalues,'ukb\wmh\wmh_corr_structural_mri/results/FA_pvalues.csv','WriteRowNames',true);



%%% ODI
ind2=strcmp(Diffusion_MRI.dti_category,'OD');
ODI_names=Diffusion_MRI(ind2,:);

ODI_rvalues=Phenotypic_dti_corr_rval(:,ind2);
ODI_pvalues=Phenotypic_dti_corr_pval(:,ind2);

ODI_names.FieldID2=ODI_rvalues.Properties.VariableNames';

ODI_rvalues.Properties.VariableNames=ODI_names.roi_names;
ODI_pvalues.Properties.VariableNames=ODI_names.roi_names;

writetable(ODI_rvalues,'ukb\wmh\wmh_corr_structural_mri/results/ODI_rvalues.csv','WriteRowNames',true);
writetable(ODI_pvalues,'ukb\wmh\wmh_corr_structural_mri/results/ODI_pvalues.csv','WriteRowNames',true);


%%% MD
ind2=strcmp(Diffusion_MRI.dti_category,'MD');
MD_names=Diffusion_MRI(ind2,:);

MD_rvalues=Phenotypic_dti_corr_rval(:,ind2);
MD_pvalues=Phenotypic_dti_corr_pval(:,ind2);

MD_names.FieldID2=MD_rvalues.Properties.VariableNames';

MD_rvalues.Properties.VariableNames=MD_names.roi_names;
MD_pvalues.Properties.VariableNames=MD_names.roi_names;

writetable(MD_rvalues,'ukb\wmh\wmh_corr_structural_mri/results/MD_rvalues.csv','WriteRowNames',true);
writetable(MD_pvalues,'ukb\wmh\wmh_corr_structural_mri/results/MD_pvalues.csv','WriteRowNames',true);


%%% ICVF
ind2=strcmp(Diffusion_MRI.dti_category,'ICVF');
ICVF_names=Diffusion_MRI(ind2,:);

ICVF_rvalues=Phenotypic_dti_corr_rval(:,ind2);
ICVF_pvalues=Phenotypic_dti_corr_pval(:,ind2);

ICVF_names.FieldID2=ICVF_rvalues.Properties.VariableNames';

ICVF_rvalues.Properties.VariableNames=ICVF_names.roi_names;
ICVF_pvalues.Properties.VariableNames=ICVF_names.roi_names;

writetable(ICVF_rvalues,'ukb\wmh\wmh_corr_structural_mri/results/ICVF_rvalues.csv','WriteRowNames',true);
writetable(ICVF_pvalues,'ukb\wmh\wmh_corr_structural_mri/results/ICVF_pvalues.csv','WriteRowNames',true);



%%% ISOVF
ind2=strcmp(Diffusion_MRI.dti_category,'ISOVF');
ISOVF_names=Diffusion_MRI(ind2,:);

ISOVF_rvalues=Phenotypic_dti_corr_rval(:,ind2);
ISOVF_pvalues=Phenotypic_dti_corr_pval(:,ind2);

ISOVF_names.FieldID2=ISOVF_rvalues.Properties.VariableNames';

ISOVF_rvalues.Properties.VariableNames=ISOVF_names.roi_names;
ISOVF_pvalues.Properties.VariableNames=ISOVF_names.roi_names;

writetable(ISOVF_rvalues,'ukb\wmh\wmh_corr_structural_mri/results/ISOVF_rvalues.csv','WriteRowNames',true);
writetable(ISOVF_pvalues,'ukb\wmh\wmh_corr_structural_mri/results/ISOVF_pvalues.csv','WriteRowNames',true);

