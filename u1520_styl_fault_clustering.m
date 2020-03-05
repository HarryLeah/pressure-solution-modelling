clear
close all

% This script produces, amongst other things, the basis for the
% frequency plots of stylolites and faults used in 'Mixed brittle and 
% ductile strain localisation and weakening in pelagic sediments seaward 
% of the Hikurangi margin, New Zealand', submitted to AGU Tectonics

% for more information view the github repository at:
% https://github.com/HarryLeah/pressure-solution-modelling
% or contact Harry Leah at LeahHR@cardiff.ac.uk

%% choose your parameters
%choose depth range for plot in m
min_depth=500;
max_depth=850;

%choose histogram bin width in m
width=1;

%threshhold bounds for variation coefficient calculation in m
thresh_min=0.0001;
thresh_max=100;

%number of threshold levels to do
thresh_n=5000;

%width of fault core in cm
core_width=5.5;

%choose distance over which to count cumulative sums in m
cum_dist=1;

%% import variables

%import samples
samples=xlsread('./data/input_samples.xlsx');

samples_top=samples(:,10);
samples_bottom=samples(:,11);
samples_avg=(samples_top+samples_bottom)./2;

formatspec="Core section %dR%dW offset %d to %d cm";
sample_labels=strings(15,1);
%labels
for i=1:15
    
    sample_labels(i)=string(sprintf(formatspec,...
        samples(i,4),...
        samples(i,6),...
        samples(i,8),...
        samples(i,9)));
    
end

%import stylolites
styl_T=xlsread('./data/stylolites_domain_2.xlsx');
styl_depth=styl_T(:,6);
styl_inten=styl_T(:,7);

%import faults
fault_T=xlsread('./data/faults_domain_2.xlsx');
fault_depth=fault_T(:,6);
fault_range=(fault_T(:,4)-fault_T(:,3))./10;

% concatenate depth arrays for joint analysis
joint_depth = cat(1,styl_depth,fault_depth);

%sort faults and stylolites
%styls
[styl_depth,styl_I]=sort(styl_depth);
styl_inten=styl_inten(styl_I);

%faults
[fault_depth,fault_I]=sort(fault_depth);
fault_range=fault_range(fault_I);

%together
joint_depth = sort(joint_depth);

%calculate apparent dips for faults
tan_fault_dip=fault_range./core_width;
fault_dip=atand(tan_fault_dip);

%import and sort calcium carbonate content
calc_carb_T=csvread('./data/depth_carb.csv');
carb_depth=calc_carb_T(:,1);
calc_carb=calc_carb_T(:,2);
[carb_depth,I]=sort(carb_depth);
calc_carb=calc_carb(I);

%calculate carbonate content line gradient

carb_grad(:,1:numel(calc_carb)-1)=nan;
carb_grad_depths(:,1:numel(calc_carb)-1)=nan;

for i=1:numel(calc_carb)-1
    
    carb_grad(i)=calc_carb(i+1)-calc_carb(i);
    
    carb_grad_depths(i)=carb_depth(i)+((carb_depth(i+1)-carb_depth(i))/2);
    
end

%import and sort MAD
MAD_T=readtable('MAD_17_10_2018_no_heds.csv');
por_depth=table2array(MAD_T(:,9));
por=table2array(MAD_T(:,17));
[por_depth,I]=sort(por_depth);
por=por(I);

%calculate stylolite spacings
total_space=max_depth-min_depth;
styl_space(1:numel(styl_depth),1)=NaN;
styl_space(1)=styl_depth(1)-min_depth;
for i=2:numel(styl_depth)
    
    if i~=numel(styl_depth)
        styl_space(i)=styl_depth(i+1)-styl_depth(i);
    else
        styl_space(i)=max_depth-styl_depth(i);
    end
    
end

%calculate fault spacings
fault_space(1:numel(fault_depth),1)=NaN;
fault_space(1)=fault_depth(1)-min_depth;
for i=2:numel(fault_depth)
    
    if i~=numel(fault_depth)
        fault_space(i)=fault_depth(i+1)-fault_depth(i);
    else
        fault_space(i)=max_depth-fault_depth(i);
    end
    
end

%calculate joint fault and stylolite spacings
joint_space(1:numel(joint_depth),1)=NaN;
joint_space(1)=joint_depth(1)-min_depth;
for i=2:numel(joint_depth)
    
    if i~=numel(joint_depth)
        joint_space(i)=joint_depth(i+1)-joint_depth(i);
    else
        joint_space(i)=max_depth-joint_depth(i);
    end
    
end



%calculate coefficient of variation for faults and styls
styl_SD=std(styl_space);
styl_s_mean=mean(styl_space);
coef_v_styl_all=styl_SD/styl_s_mean;

fault_SD=std(fault_space);
fault_s_mean=mean(fault_space);
coef_v_fault_all=fault_SD/fault_s_mean;

joint_SD=std(joint_space);
joint_s_mean=mean(joint_space);
coef_v_joint_all=joint_SD/joint_s_mean;

%create plottable line
fault_cv_all_plot(1)=coef_v_fault_all;
fault_cv_all_plot(2)=coef_v_fault_all;

styl_cv_all_plot(1)=coef_v_styl_all;
styl_cv_all_plot(2)=coef_v_styl_all;

joint_cv_all_plot(1)=coef_v_joint_all;
joint_cv_all_plot(2)=coef_v_joint_all;

cv_plot_bounds(1)=thresh_min;
cv_plot_bounds(2)=thresh_max;

%thresholding and recalculation of Cv
thresh_min=log10(thresh_min);
thresh_max=log10(thresh_max);
thresh=logspace(thresh_min,thresh_max,thresh_n);

coef_v_styl(1:numel(thresh))=NaN;
coef_v_fault(1:numel(thresh))=NaN;
coef_v_joint(1:numel(thresh))=NaN;

%coef_v_styl(1)=coef_v_styl_all;
%coef_v_fault(1)=coef_v_fault_all;

styl_space_thresh=styl_space;
fault_space_thresh=fault_space;
joint_space_thresh=joint_space;

for i=1:numel(thresh)
    
    styl_index=(styl_space_thresh>=thresh(i)); 
    fault_index=(fault_space_thresh>=thresh(i)); 
    joint_index=(joint_space_thresh>=thresh(i));
    
    styl_space_thresh=styl_space_thresh(styl_index);
    fault_space_thresh=fault_space_thresh(fault_index);
    joint_space_thresh=joint_space_thresh(joint_index);
    
    styl_SD=std(styl_space_thresh);
    styl_s_mean=mean(styl_space_thresh);
    coef_v_styl(i)=styl_SD/styl_s_mean;
    
    fault_SD=std(fault_space_thresh);
    fault_s_mean=mean(fault_space_thresh);
    coef_v_fault(i)=fault_SD/fault_s_mean;
    
    joint_SD=std(fault_space_thresh);
    joint_s_mean=mean(joint_space_thresh);
    coef_v_joint(i)=joint_SD/joint_s_mean;
    
    
end

%calculate histcounts
edges=min_depth:width:max_depth;
s_counts=histcounts(styl_depth,edges);
f_counts=histcounts(fault_depth,edges);
f_bin_locs=(min_depth+(width/2)):width:(max_depth);
s_bin_locs=(min_depth+(width/2)):width:(max_depth);

ylimmm=([min_depth max_depth]);

%% calculate Cv for lithologies
major_litho=xlsread('./data/Simplified_domain_2_strat.xlsx');

strat_size=size(major_litho);

litho_thickness=major_litho(:,2)-major_litho(:,1);

%loop should be fine as long as thresh, coef_v_styl, and coef_v_fault, cv_plot_bounds, and styl_cv_all_plot aren't touched

for z=1:strat_size(1)
    
    litho_top=major_litho(z,1);
    litho_base=major_litho(z,2);
    
    litho_styl_index=(styl_depth>=litho_top) & (styl_depth<=litho_base);
    litho_fault_index=(fault_depth>=litho_top) & (fault_depth<=litho_base);
    litho_joint_index=(joint_depth>=litho_top) & (joint_depth<=litho_base);
    
    litho_styl_depth=styl_depth(litho_styl_index);
    litho_fault_depth=fault_depth(litho_fault_index);
    litho_joint_depth=joint_depth(litho_joint_index);
    
    styl_numbers(z)=numel(litho_styl_depth);
    fault_numbers(z)=numel(litho_fault_depth);
    joint_numbers(z)=numel(litho_joint_depth);
    
    if numel(litho_styl_depth)~=0
    if numel(litho_fault_depth)~=0
    if numel(litho_joint_depth)~=0
    
    %calculate stylolite spacings
    litho_space=litho_base-litho_top;
    litho_styl_space(1:numel(litho_styl_depth),1)=NaN;
    litho_styl_space(1)=litho_styl_depth(1)-litho_top;
    
for i=2:numel(litho_styl_depth)
    
    if i~=numel(litho_styl_depth)
        litho_styl_space(i)=litho_styl_depth(i+1)-litho_styl_depth(i);
    else
        litho_styl_space(i)=litho_base-litho_styl_depth(i);
    end
    
end

    %calculate fault spacings
    litho_fault_space(1:numel(litho_fault_depth),1)=NaN;
    litho_fault_space(1)=litho_fault_depth(1)-litho_top;

for i=2:numel(litho_fault_depth)
    
    if i~=numel(litho_fault_depth)
        litho_fault_space(i)=litho_fault_depth(i+1)-litho_fault_depth(i);
    else
        litho_fault_space(i)=litho_base-litho_fault_depth(i);
    end
    
end

   %calculate joint fault and stylolite spacings
    litho_joint_space(1:numel(litho_joint_depth),1)=NaN;
    litho_joint_space(1)=litho_joint_depth(1)-litho_top;

for i=2:numel(litho_joint_depth)
    
    if i~=numel(litho_joint_depth)
        litho_joint_space(i)=litho_joint_depth(i+1)-litho_joint_depth(i);
    else
        litho_joint_space(i)=litho_base-litho_joint_depth(i);
    end
    
end

litho_styl_space_thresh=litho_styl_space;
litho_fault_space_thresh=litho_fault_space;
litho_joint_space_thresh=litho_joint_space;

for i=1:numel(thresh)
    
    litho_styl_index=(litho_styl_space_thresh>=thresh(i)); 
    litho_fault_index=(litho_fault_space_thresh>=thresh(i)); 
    litho_joint_index=(litho_joint_space_thresh>=thresh(i)); 
    
    litho_styl_space_thresh=litho_styl_space_thresh(litho_styl_index);
    litho_fault_space_thresh=litho_fault_space_thresh(litho_fault_index);
    litho_joint_space_thresh=litho_joint_space_thresh(litho_joint_index);
    
    litho_styl_SD=std(litho_styl_space_thresh);
    litho_styl_s_mean=mean(litho_styl_space_thresh);
    litho_coef_v_styl(i,z)=litho_styl_SD/litho_styl_s_mean;
    
    litho_fault_SD=std(litho_fault_space_thresh);
    litho_fault_s_mean=mean(litho_fault_space_thresh);
    litho_coef_v_fault(i,z)=litho_fault_SD/litho_fault_s_mean;
    
    litho_joint_SD=std(litho_joint_space_thresh);
    litho_joint_s_mean=mean(litho_joint_space_thresh);
    litho_coef_v_joint(i,z)=litho_joint_SD/litho_joint_s_mean;
    
end

    end
    end
    end
end

%calculate number of features per m on average
styl_counts_p_m=rot90(styl_numbers,3)./litho_thickness;
fault_counts_p_m=rot90(fault_numbers,3)./litho_thickness;
joint_counts_p_m=rot90(joint_numbers,3)./litho_thickness;

%% calculate cumulative counts and gradient changes

%for faults and stylolites
cum_fault(:,1:ceil(((max_depth-min_depth)/cum_dist)))=nan;
cum_styl(:,1:ceil(((max_depth-min_depth)/cum_dist)))=nan;
cum_joint(:,1:ceil(((max_depth-min_depth)/cum_dist)))=nan;

for i=1:ceil(((max_depth-min_depth)/cum_dist)+1)
    
    cum_fault(i)=sum(fault_depth(:)<=(min_depth+(i*cum_dist))); 
    
    cum_styl(i)=sum(styl_depth(:)<=(min_depth+(i*cum_dist))); 

    cum_joint(i)=sum(joint_depth(:)<=(min_depth+(i*cum_dist))); 
    
end

%calculate gradients
cum_fault_grad(:,1:numel(cum_fault)-1)=nan;
cum_styl_grad(:,1:numel(cum_styl)-1)=nan;
cum_joint_grad(:,1:numel(cum_joint)-1)=nan;
cum_grad_depth(:,1:numel(cum_styl)-1)=nan;
depths=min_depth:cum_dist:max_depth;

for i=1:numel(cum_fault)-1
    
    cum_fault_grad(i)=cum_fault(i+1)-cum_fault(i);
    cum_styl_grad(i)=cum_styl(i+1)-cum_styl(i);
    cum_joint_grad(i)=cum_joint(i+1)-cum_joint(i);
    
    cum_grad_depth(i)=depths(i)+(cum_dist/2)+cum_dist;
    
end

%Normalising threshold thicknesses to unit thicknesses

thresh_calc_mud=thresh./litho_thickness(1);
thresh_marl_upp=thresh./litho_thickness(2);
thresh_mud=thresh./litho_thickness(3);
thresh_marl_low=thresh./litho_thickness(4);
thresh_chalk_upp=thresh./litho_thickness(5);
thresh_debris=thresh./litho_thickness(6);
thresh_chalk_low=thresh./litho_thickness(7);
thresh_all_norm=thresh./(max_depth-min_depth);

thresh_litho_comb(:,1)=thresh_calc_mud;
thresh_litho_comb(:,2)=thresh_marl_upp;
thresh_litho_comb(:,3)=thresh_mud;
thresh_litho_comb(:,4)=thresh_marl_low;
thresh_litho_comb(:,5)=thresh_chalk_upp;
thresh_litho_comb(:,6)=thresh_debris;
thresh_litho_comb(:,7)=thresh_chalk_low;
thresh_litho_comb(:,8)=thresh_all_norm;


for i=1:8
    
    if i==8
        one_m_x(i)=1/(max_depth-min_depth);
        
        one_m_y_styl(i)=interp1(thresh_litho_comb(:,i),coef_v_styl,one_m_x(i));
        
        one_m_y_fault(i)=interp1(thresh_litho_comb(:,i),coef_v_fault,one_m_x(i));
        
        one_m_y_joint(i)=interp1(thresh_litho_comb(:,i),coef_v_joint,one_m_x(i));
        
    else
        one_m_x(i)=1/litho_thickness(i);
            
        %interpolate value for stylolite
        one_m_y_styl(i)=interp1(thresh_litho_comb(:,i),litho_coef_v_styl(:,i),one_m_x(i));
        
        %interpolate value for faults
        one_m_y_fault(i)=interp1(thresh_litho_comb(:,i),litho_coef_v_fault(:,i),one_m_x(i));
        
        %interpolate value for joint faults and stylolites
        one_m_y_joint(i)=interp1(thresh_litho_comb(:,i),litho_coef_v_joint(:,i),one_m_x(i));
    end

    
end



%% plotting coefficient of variation

figure
title('Coefficient of variation feature-space length-scaling for domain 2')
plot(thresh,coef_v_styl,'-k')
ylabel('Coefficient of variation')
xlabel('Threshold thickness (m)')
set(gca,'XScale','log')
xlim([min(thresh) max(thresh)])
hold on
plot(thresh,coef_v_fault,'-r')
plot(thresh,coef_v_joint,'-b')

plot(cv_plot_bounds,styl_cv_all_plot,'--k')
plot(cv_plot_bounds,fault_cv_all_plot,'--r')
plot(cv_plot_bounds,joint_cv_all_plot,'--b')

legend('Stylolite C_v','Fault C_v','Stylolite and fault C_v','Entire sample stylolite C_v','Entire sample fault C_v','Entire sample stylolite and fault C_v')

hold off


%% plotting cumulative plot with gradient

%scatter of faults and stylolites
figure
title('Frequency of stylolites and faults down borehole')
subplot(1,4,1)
scatter(styl_inten,styl_depth,20,'k','filled')
ylabel('Depth (mbsf)')
xlabel('Stylolite intensity')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([0 5])
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
scatter(fault_dip,fault_depth,20,'r','filled')
xlabel('Fault apparent dip (degrees)')
set(gca,'YDir','reverse')
set(gca,'YColor','none')
ylim(ylimmm)
xlim([0 90])
set(gca,'XTick',0:30:90)
hold off

%plot of cumulative frequency and gradient for stylolites
subplot(1,4,2)
plot(cum_styl,depths,'-b')
ylabel('Depth (mbsf)')
xlabel('Cumulative stylolite number')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([-100 1000])
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
plot(cum_styl_grad,cum_grad_depth,'-k')
xlabel('Gradient of cumulative stylolite number')
set(gca,'YDir','reverse')
set(gca,'YColor','none')
ylim(ylimmm)
xlim([-5 80])
hold off

%plot of cumulative frequency and gradient for faults
subplot(1,4,3)
plot(cum_fault,depths,'-r')
ylabel('Depth (mbsf)')
xlabel('Cumulative fault number')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([-40 400])
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
plot(cum_fault_grad,cum_grad_depth,'-k')
xlabel('Gradient of cumulative fault number')
set(gca,'YDir','reverse')
set(gca,'YColor','none')
ylim(ylimmm)
xlim([-1 20])
hold off

%plot CaCO3 content and gradient
subplot(1,4,4)
plot(calc_carb,carb_depth,'-k')
xlabel('CaC0_3 Weight percent')
set(gca,'YDir','reverse')
ylim(ylimmm)

figure
%plot of cumulative frequency and gradient for joint
plot(cum_joint,depths,'-r')
ylabel('Depth (mbsf)')
xlabel('Cumulative joint number')
set(gca,'YDir','reverse')
ylim(ylimmm)
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
plot(cum_joint_grad,cum_grad_depth,'-k')
xlabel('Gradient of cumulative joint number')
set(gca,'YDir','reverse')
set(gca,'YColor','none')
ylim(ylimmm)
hold off
%% plotting Cv by lithology

figure
plot(thresh,litho_coef_v_fault(:,1),...
    thresh,litho_coef_v_fault(:,2),...
    thresh,litho_coef_v_fault(:,3),...
    thresh,litho_coef_v_fault(:,4),...
    thresh,litho_coef_v_fault(:,5),...
    thresh,litho_coef_v_fault(:,6),...
    thresh,litho_coef_v_fault(:,7))
ylabel('Coefficient of variation')
xlabel('Threshold thickness (m)')
title('Coefficient of variation fault-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh) max(thresh)])
ylim([-0.2 4])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Mudstone 720.93-738.68	mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf')

figure
plot(thresh,litho_coef_v_joint(:,1),...
    thresh,litho_coef_v_joint(:,2),...
    thresh,litho_coef_v_joint(:,3),...
    thresh,litho_coef_v_joint(:,4),...
    thresh,litho_coef_v_joint(:,5),...
    thresh,litho_coef_v_joint(:,6),...
    thresh,litho_coef_v_joint(:,7))
ylabel('Coefficient of variation')
xlabel('Threshold thickness (m)')
title('Coefficient of variation fault and stylolite-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh) max(thresh)])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Mudstone 720.93-738.68	mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf')


figure
plot(thresh,litho_coef_v_styl(:,1),...
    thresh,litho_coef_v_styl(:,2),...
    thresh,litho_coef_v_styl(:,4),...
    thresh,litho_coef_v_styl(:,5),...
    thresh,litho_coef_v_styl(:,6),...
    thresh,litho_coef_v_styl(:,7))
ylabel('Coefficient of variation')
xlabel('Threshold thickness (m)')
title('Coefficient of variation stylolite-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh) max(thresh)])
ylim([-0.2 6.5])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf')

%% bar charts of number per metre by lithology
lithos=categorical(...
{
    '510.96-655.15 mbsf Calcareous mudstone',...
    '655.15-720.93 mbsf Marl',...
    '720.93-738.68 mbsf Mudstone',...
    '738.68-773.26 mbsf Marl',...
    '773.26-781.12 mbsf Chalk',...
    '781.12-797.99 mbsf Debris flow',...
    '797.99-848.45 mbsf Chalk'
});

figure
subplot(1,2,1)
barh(lithos,styl_counts_p_m,'FaceColor','blue')
set(gca,'YDir','reverse')
xlabel('Stylolite counts per metre')

subplot(1,2,2)
barh(lithos,fault_counts_p_m,'FaceColor','red')
set(gca,'YDir','reverse')
xlabel('Fault counts per metre')
xlim([0 4])

%% plot combined stylolite diagram
figure
subplot(1,2,1)
plot(cum_styl,depths,'-b')
ylabel('Depth (mbsf)')
xlabel('Cumulative stylolite number')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([-99 1000])
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
plot(calc_carb,carb_depth,'-k',cum_styl_grad,cum_grad_depth,'-b')
xlabel('CaC0_3 Weight percent (black) or gradient of stylolite cumulative frequency (blue)')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([0 100])
hold off

subplot(1,2,2)
scatter(styl_inten,styl_depth,20,'k','filled')
ylabel('Depth (mbsf)')
xlabel('Stylolite intensity')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([0 5])


%% plot combined fault diagram

figure
subplot(1,2,1)
plot(cum_fault,depths,'-r')
ylabel('Depth (mbsf)')
xlabel('Cumulative fault number (red)')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([-49 400])

subplot(1,2,2)
plot(cum_fault_grad,cum_grad_depth,'-k')
xlabel('Gradient of Cumulative fault number (black)')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([0 20])
ax1 = gca;
box(ax1,'off')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none');
hold on
scatter(fault_dip,fault_depth,20,'r','filled')
ylabel('Depth (mbsf)')
xlabel('Apparent fault dip')
set(gca,'YDir','reverse')
ylim(ylimmm)
xlim([0 90])
hold off

%% plot normalised Cv plots by lithology

figure
plot(thresh_calc_mud,litho_coef_v_fault(:,1),...
    thresh_marl_upp,litho_coef_v_fault(:,2),...
    thresh_mud,litho_coef_v_fault(:,3),...
    thresh_marl_low,litho_coef_v_fault(:,4),...
    thresh_chalk_upp,litho_coef_v_fault(:,5),...
    thresh_debris,litho_coef_v_fault(:,6),...
    thresh_chalk_low,litho_coef_v_fault(:,7),...
    thresh_all_norm,coef_v_fault,'-k')
ylabel('Coefficient of variation')
xlabel('Normalised threshold thickness (threshold thickness/unit thickness)')
title('Normalised coefficient of variation fault-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh_calc_mud) 1])
ylim([-0.2 4.7])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Mudstone 720.93-738.68	mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf',...
    'All faults')
hold on
for i=1:8
    if i==8
        scatter(one_m_x(i),one_m_y_fault(i),20,...
            'k','filled','DisplayName','1 m spacing');
    else
        scatter(one_m_x(i),one_m_y_fault(i),20,...
            'filled','DisplayName','1 m spacing');
    end
end
hold off

figure
plot(thresh_calc_mud,litho_coef_v_styl(:,1),...
    thresh_marl_upp,litho_coef_v_styl(:,2),...
    thresh_marl_low,litho_coef_v_styl(:,4),...
    thresh_chalk_upp,litho_coef_v_styl(:,5),...
    thresh_debris,litho_coef_v_styl(:,6),...
    thresh_chalk_low,litho_coef_v_styl(:,7),...
    thresh_all_norm,coef_v_styl,'-k')
ylabel('Coefficient of variation')
xlabel('Normalised threshold thickness (threshold thickness/unit thickness)')
title('Normalised coefficient of variation stylolite-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh_calc_mud) 1])
ylim([-0.2 7])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf',...
    'All stylolites')
hold on
for i=1:8
    if i==8
        scatter(one_m_x(i),one_m_y_styl(i),20,...
            'k','filled','DisplayName','1 m spacing');
    elseif i==3
        
    else
        scatter(one_m_x(i),one_m_y_styl(i),20,...
            'filled','DisplayName','1 m spacing');
    end
end
hold off

figure
map = viridisColorMap(7);
threshes = [thresh_calc_mud;thresh_marl_upp;thresh_mud;thresh_marl_low;thresh_chalk_upp;thresh_debris;thresh_chalk_low];
hold on
for i=1:7
plot(threshes(i,:),litho_coef_v_joint(:,i),'Color',map(i,:),'LineWidth',1.5)

end
plot(thresh_all_norm,coef_v_joint,'-k','LineWidth',1.5)
ylabel('Coefficient of variation')
xlabel('Normalised threshold thickness (threshold thickness/unit thickness)')
title('Normalised coefficient of variation in joint stylolite and fault-gap length-scaling by lithology')
set(gca,'XScale','log')
xlim([min(thresh_calc_mud) 1])
legend('Calcareous mudstone 510.96-655.15 mbsf',...
    'Marl 655.15-720.93 mbsf',...
    'Mudstone 720.93-738.68	mbsf',...
    'Marl 738.68-773.26 mbsf',...
    'Chalk 773.26-781.12 mbsf',...
    'Debris flow 781.12-797.99 mbsf',...
    'Chalk 797.99-848.45 mbsf',...
    'All faults and stylolites')
for i=1:8
    if i==8
        scatter(one_m_x(i),one_m_y_joint(i),20,...
            'filled','DisplayName','1 m spacing',...
            'MarkerFaceColor','k',...
            'MarkerEdgeColor','k');
    else
        scatter(one_m_x(i),one_m_y_joint(i),20,...
            'filled','DisplayName','1 m spacing',...
            'MarkerFaceColor',map(i,:),...
            'MarkerEdgeColor',map(i,:));
    end
end
hold off


%% compare fault and stylolite clustering

figure
plot(cum_fault_grad.*4,cum_grad_depth,'-r',cum_styl_grad,cum_grad_depth,'-b')
xlabel('Gradient of Cumulative fault number (black)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse')
ylim(ylimmm)
legend('Gradient of cumulative fault number*4','Gradient of cumulative stylolite number','Location','northeast')

%% plot sample depths
clear litho_bounds
%create double value for lines
litho_bounds(:,1)=major_litho(:,1);
litho_bounds(:,2)=major_litho(:,1);
litho_bounds(8,1)=major_litho(7,2);
litho_bounds(8,2)=major_litho(7,2);

x_axis=[0 100];

figure
plot(cum_fault_grad.*4,cum_grad_depth,'-r',cum_styl_grad,cum_grad_depth,'-b')
xlabel('Gradient of Cumulative fault number (black)')
ylabel('Depth (mbsf)')
hold on
for i=1:8
    plot(x_axis,litho_bounds(i,:),'-k')
    
    for n=1:7
        text(50,((litho_bounds(n,1)+litho_bounds(n+1,1))/2),char(lithos(n)))
    end
    
    
end
ylim([650 max_depth])
set(gca,'YDir','reverse')

for i=1:15
    
    scatter(0.45*100,samples_avg(i),'filled')
    text(0.5*100,samples_avg(i),sample_labels(i))
    
end

hold off

%% plot histograms
%plot dip histogram
figure
dip_plot=polaraxes;
polarhistogram(deg2rad(fault_dip),deg2rad([0:5:90]),'FaceColor','r')
dip_plot.ThetaZeroLocation='right';
dip_plot.ThetaDir='clockwise';
thetalim([0 90])
dip_plot.FontSize=15;


%plot stylolite intensity histogram
figure
histogram(styl_inten,'FaceColor','b')
ax=gca;
ax.YScale='log';
ax.FontSize=15;
yticks([1e0 1e1 1e2 1e3])
yticklabels([1 10 100 1000])
xlim([0.5 5.5])
xticks([1 2 3 4 5])