
close all
clear all
clc



cd  ukb\wmh\wmh_corr_structural_mri

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
Diffusion_MRI.side=side;


% %
ukb_wmh=readtable('ukb\wmh\YJ/s5_Pred.csv');
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

[final_eid,~]= intersect(intersect(intersect(intersect(ukb_wmh.eid,Diffusion_MRI_bl_data.eid),wmh_cov_info_0.eid),wmh_sur_cov.eid),cov_info_2.eid);



%%
[~,ind2]= intersect(Diffusion_MRI_bl_data.eid,final_eid);
Diffusion_MRI_bl_data1=Diffusion_MRI_bl_data(ind2,:);


[~,ind2]= intersect(ukb_wmh.eid,final_eid);
ukb_wmh=ukb_wmh(ind2,:);



[~,ind2]= intersect(wmh_cov_info_0.eid,final_eid);
wmh_cov_info_0=wmh_cov_info_0(ind2,:);

[~,ind2]= intersect(wmh_sur_cov.eid,final_eid);
wmh_sur_cov=wmh_sur_cov(ind2,:);

[~,ind2]= intersect(cov_info_2.eid,final_eid);
cov_info_2=cov_info_2(ind2,:);


data_eid=[ukb_wmh.eid,Diffusion_MRI_bl_data1.eid,wmh_cov_info_0.eid];

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

%instance2_duration=zscore(instance2_duration./12);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


instance2_age=wmh_cov_info_0.age+instance2_duration./12;


% 26518	'Volume of TotalGray (whole brain)'
instance2_cov_info=[wmh_cov_info_0.sex,zscore(instance2_age), bl_Ethnic,instance2_site,bl_Medical_covinfo];%,
%[T,Res,total_wmh_res,beta]=BWAS_Tregression(instance2_cov_info,wmh_log2);



wmh_corr_dkt_results=nan(size(Diffusion_MRI_names,1),2);
corr_rval=nan(size(Diffusion_MRI_names,1),1);
corr_pval=nan(size(Diffusion_MRI_names,1),1);

ValueType=[];
for i=1:size(Diffusion_MRI_names,1)

    tic

    x1=Diffusion_MRI_bl_data1(:,i+1).Variables;
    ind2=~isnan(x1);

    x=x1(ind2,:);
    y=ukb_wmh(ind2,:).y_pred;
    z=instance2_cov_info(ind2,:);

    [corr_rval(i,1),corr_pval(i,1)]= partialcorr(x,y,z);% );%

    toc

end



wmh_corr_dkt_results=array2table([corr_rval,corr_pval],"VariableNames",{'rvalue','raw_p'},'RowNames', Diffusion_MRI_names);
wmh_corr_dkt_results.Diffusion_MRI_names=Diffusion_MRI_names;

%wmh_corr_dkt_results.BF_P=corr_pval .* length(corr_pval);


Diffusion_MRI.rvalue=wmh_corr_dkt_results.rvalue;
Diffusion_MRI.pvalue=wmh_corr_dkt_results.raw_p;


ind2=strcmp(Diffusion_MRI.dti_category,'FA');
FA_Diffusion_MRI=Diffusion_MRI(ind2,:);

ind2=strcmp(Diffusion_MRI.dti_category,'MD');
MD_Diffusion_MRI=Diffusion_MRI(ind2,:);

ind2=strcmp(Diffusion_MRI.dti_category,'ISOVF');
ISOVF_Diffusion_MRI=Diffusion_MRI(ind2,:);

ind2=strcmp(Diffusion_MRI.dti_category,'ICVF');
ICVF_Diffusion_MRI=Diffusion_MRI(ind2,:);

ind2=strcmp(Diffusion_MRI.dti_category,'OD');
OD_Diffusion_MRI=Diffusion_MRI(ind2,:);


%%
roi_names_dti=[FA_Diffusion_MRI.roi_names,MD_Diffusion_MRI.roi_names,ISOVF_Diffusion_MRI.roi_names,ICVF_Diffusion_MRI.roi_names];


pre_wmh_dti_r_array=[FA_Diffusion_MRI.rvalue,MD_Diffusion_MRI.rvalue,OD_Diffusion_MRI.rvalue ,ICVF_Diffusion_MRI.rvalue,ISOVF_Diffusion_MRI.rvalue];
pre_wmh_dti_r_table=array2table(pre_wmh_dti_r_array,"RowNames",FA_Diffusion_MRI.roi_names,"VariableNames",{'FA','MD','ODI','ICVF','ISOVF'});

pre_wmh_dti_p_array=[FA_Diffusion_MRI.pvalue,MD_Diffusion_MRI.pvalue,OD_Diffusion_MRI.pvalue ,ICVF_Diffusion_MRI.pvalue,ISOVF_Diffusion_MRI.pvalue];
pre_wmh_dti_p_array=pre_wmh_dti_p_array*(27*5);
pre_wmh_dti_p_table=array2table(pre_wmh_dti_p_array,"RowNames",FA_Diffusion_MRI.roi_names,"VariableNames",{'FA','MD','ODI','ICVF','ISOVF'});

writetable(pre_wmh_dti_r_table,['ukb/wmh/wmh_corr_structural_mri/results/pre_wmh_dti_r_table.csv'],'WriteRowNames',true);
writetable(pre_wmh_dti_p_table,['ukb/wmh/wmh_corr_structural_mri/results/pre_wmh_dti_p_table.csv'],'WriteRowNames',true);


%%
pre_wmh_dti_r_table=readtable('ukb/wmh/wmh_corr_structural_mri/results/pre_wmh_dti_r_table.csv');



