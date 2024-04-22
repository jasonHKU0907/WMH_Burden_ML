close all
clear all
clc

addpath(genpath('software\ENIGMA-2.0.0'))
load ukb\wmh\ukb_subsets/Freesurfer_DKT_data.mat


% %
structural_type={'Mean intensity of ','Mean thickness of ','Volume of ','Area of '};


Freesurfer_DKT.Description2=Freesurfer_DKT.Description;
Freesurfer_DKT.type_ind2=zeros(size(Freesurfer_DKT,1),1);

for k=1:size(Freesurfer_DKT,1)
    t1_name= Freesurfer_DKT(k,:).Description;

    for n=1:length(structural_type)

        if contains(Freesurfer_DKT(k,:).Description,structural_type{n})
            Freesurfer_DKT(k,:).type_ind2=n;
            t2_name=strrep(t1_name,structural_type{n},'');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t2_name=lower(t2_name);

    if contains(t2_name,'right')
        t3_name=strrep(t2_name,' (right hemisphere)','');
        t3_name=strcat('R_',t3_name);

    elseif contains(t2_name,'left')
        t3_name=strrep(t2_name,' (left hemisphere)','');
        t3_name=strcat('L_',t3_name);
    else
        t3_name=t2_name;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Freesurfer_DKT(k,:).Description2=t3_name;
end

% sum(Freesurfer_DKT.type_ind2==2)
% 



wmh_corr_dkt_results=readtable('ukb/wmh/wmh_corr_structural_mri/results/predict_wmh_corr_dkt_results.csv');
wmh_corr_dkt_results.Description2=Freesurfer_DKT.Description2;

wmh_corr_dkt_results.adjusted_p=wmh_corr_dkt_results.raw_p * (62*3+16);

ind2=Freesurfer_DKT.type_ind2==3
wmh_corr_volume_results=wmh_corr_dkt_results(ind2,:);

ind2=Freesurfer_DKT.type_ind2==2
wmh_corr_area_results=wmh_corr_dkt_results(ind2,:);


ind2=Freesurfer_DKT.type_ind2==4
wmh_corr_thickness_results=wmh_corr_dkt_results(ind2,:);



CortThick=readtable("ENIGMA-2.0.0\matlab\shared\data\summary_statistics/mddadult_case-controls_CortThick.csv");
Cortical_order=CortThick.Structure;


%%    volume   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wmh_corr_cortical_volume=CortThick(:,2);
wmh_corr_cortical_volume.rvalues=zeros(size(Cortical_order));
wmh_corr_cortical_volume.adjusted_p=ones(size(Cortical_order));

for k=1:length(Cortical_order)
    ind2=strcmp(wmh_corr_volume_results.Description2,Cortical_order{k});
    if sum(ind2)>0
        wmh_corr_cortical_volume(k,:).rvalues= wmh_corr_volume_results(ind2,:).rvalue;
        wmh_corr_cortical_volume(k,:).adjusted_p= wmh_corr_volume_results(ind2,:).adjusted_p;
    else
        Cortical_order{k};
    end
end

cortical_volume_r=wmh_corr_cortical_volume.rvalues;
cortical_volume_r(wmh_corr_cortical_volume.adjusted_p>0.05)=0;

 
% %    area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wmh_corr_cortical_area=CortThick(:,2);
wmh_corr_cortical_area.rvalues=zeros(size(Cortical_order));
wmh_corr_cortical_area.adjusted_p=ones(size(Cortical_order));

for k=1:length(Cortical_order)
    ind2=strcmp(wmh_corr_area_results.Description2,Cortical_order{k});
    if sum(ind2)>0
        wmh_corr_cortical_area(k,:).rvalues= wmh_corr_area_results(ind2,:).rvalue;
        wmh_corr_cortical_area(k,:).adjusted_p= wmh_corr_area_results(ind2,:).adjusted_p;
    else
        Cortical_order{k};
    end
end

cortical_area_r=wmh_corr_cortical_area.rvalues;
cortical_area_r(wmh_corr_cortical_area.adjusted_p>0.05)=0;

% %    thickness   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wmh_corr_cortical_thickness=CortThick(:,2);
wmh_corr_cortical_thickness.rvalues=zeros(size(Cortical_order));
wmh_corr_cortical_thickness.adjusted_p=ones(size(Cortical_order));

for k=1:length(Cortical_order)
    ind2=strcmp(wmh_corr_thickness_results.Description2,Cortical_order{k});
    if sum(ind2)>0
        wmh_corr_cortical_thickness(k,:).rvalues= wmh_corr_thickness_results(ind2,:).rvalue;
        wmh_corr_cortical_thickness(k,:).adjusted_p= wmh_corr_thickness_results(ind2,:).adjusted_p;
    else
        Cortical_order{k};
    end
end

cortical_thickness_r=wmh_corr_cortical_thickness.rvalues;
cortical_thickness_r(wmh_corr_cortical_thickness.adjusted_p>0.05)=0;




%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% left-accumbens, left-amygdala, left-caudate, left-hippocampus, left-pallidum, left-putamen, left-thalamus, left-ventricles,
% right-accumbens, right-amygdala, right-caudate, right-hippocampus, right-pallidum, right-putamen, right-thalamus, right-ventricles};
SubVol=readtable("D:\software\ENIGMA-2.0.0\matlab\shared\data\summary_statistics/mddearly_case-controls_SubVol.csv");
Subcortical_order=SubVol.Structure;


Subcortical_order={'L_accumbens-area';'L_amygdala';'L_caudate';'L_hippocampus';'L_pallidum';'L_putamen';'L_thalamus-proper';'L_ventraldc';....
    'R_accumbens-area';'R_amygdala';'R_caudate';'R_hippocampus';'R_pallidum';'R_putamen';'R_thalamus-proper';'R_ventraldc'};%R_ventraldc  R_ventraldc



wmh_corr_subcortical_volume=[];
for k=1:length(Subcortical_order)

    ind2=strcmp(wmh_corr_volume_results.Description2,Subcortical_order{k});
    if sum(ind2)>0
        wmh_corr_subcortical_volume=[wmh_corr_subcortical_volume;wmh_corr_volume_results(ind2,:)];
    else
        Subcortical_order{k}
    end
end

wmh_corr_subcortical_volume.Subcortical_order=Subcortical_order;

%wmh_corr_subcortical_volume.adjusted_p=wmh_corr_subcortical_volume.raw_p*84;

subcortical_r=wmh_corr_subcortical_volume.rvalue;
subcortical_r(wmh_corr_subcortical_volume.adjusted_p>0.05)=0;


%%

% r_max=max([cortical_thickness_r;cortical_area_r;cortical_volume_r;subcortical_r]);
% r_min=min([cortical_thickness_r;cortical_area_r;cortical_volume_r;subcortical_r]);

r_max=max(abs([cortical_thickness_r;cortical_area_r;cortical_volume_r;subcortical_r]));

r_min= -r_max;

%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load summary statistics
% Map parcellated data to the surface
CT_d=cortical_volume_r
CT_d_fsa5 = parcel_to_surface(CT_d, 'aparc_fsa5');

% Project the results on the surface brain
f = figure(1), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'position',[500 500 800 400])
plot_cortical(CT_d_fsa5, 'surface_name', 'fsa5', 'color_range', ...
    [r_min r_max], 'cmap', 'RdBu_r')


print(figure(1), '-dtiff','-r900', ['figs/pre_wmh_corr_cortical_volume.tiff'])

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load summary statistics
% Map parcellated data to the surface
CT_d=cortical_thickness_r
CT_d_fsa5 = parcel_to_surface(CT_d, 'aparc_fsa5');

% Project the results on the surface brain
f = figure(2), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'position',[500 500 800 400])
plot_cortical(CT_d_fsa5, 'surface_name', 'fsa5', 'color_range', ...
    [r_min r_max], 'cmap', 'RdBu_r')


print(figure(2), '-dtiff','-r900', ['figs/pre_wmh_corr_cortical_thickness.tiff'])

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load summary statistics
% Map parcellated data to the surface
CT_d=cortical_area_r
CT_d_fsa5 = parcel_to_surface(CT_d, 'aparc_fsa5');

% Project the results on the surface brain
f = figure(3), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'position',[500 500 800 400])
plot_cortical(CT_d_fsa5, 'surface_name', 'fsa5', 'color_range', ...
    [r_min r_max], 'cmap', 'RdBu_r')


print(figure(3), '-dtiff','-r900', ['figs/pre_wmh_corr_cortical_area.tiff'])

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



SV_d=subcortical_r
% Project the results on the surface brain

f = figure(4),  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'position',[500 500 800 400])

plot_subcortical(SV_d, 'color_range', [r_min r_max], 'cmap', 'RdBu_r')

print(figure(4), '-dtiff','-r900', ['figs/pre_wmh_corr_subcortical_volume.tiff'])
