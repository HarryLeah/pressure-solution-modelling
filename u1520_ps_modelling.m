clear
close all

% This script produces, amongst other things, the basis for the
% frequency plots of stylolites and faults used in 'Mixed brittle and 
% ductile strain localisation and weakening in pelagic sediments seaward 
% of the Hikurangi margin, New Zealand', submitted to AGU Tectonics

% This script was produced by me during the first couple of years of my
% PhD. It is by no means the most efficient or well-written of code, but it
% should produce all the plots (and more) found in the corresponding paper.

% for more information view the github repository at:
% https://github.com/HarryLeah/pressure-solution-modelling
% or contact Harry Leah at LeahHR@cardiff.ac.uk

%% Set figure parameters
set(groot,'defaultFigurePaperSize',[29.7000 21])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
map_style = 'imola';

%% read Moisture and Density value tables
%load values from excel spreadsheet
S1 = xlsread('./data/MAD_12_2_2019.xlsx');

fd_fontsize = 17;


%define values from table
depth=S1(:,9);      %depth in m
bulk=S1(:,14);       %bulk density in g/cm3
water_depth=3522.1; %water depth in m
water_dens=1.024;   %water density in g/cm3
g=9.80665;          %grav acceleration in m/s
porosity=S1(:,17);

% % conversions
bulk=bulk.*1e3;     %convert bulk to kg m3
water_dens=water_dens.*1e3; %convert water to kg m3

%calculate stress sigma=Pgz
%P=density
%g=gravitational constant
%z=depth
%stress due to water depth

sigma_w=water_depth*water_dens*g;

hydrostatic=(depth.*water_dens.*g)+sigma_w;

sigma_bulk_no_w(1:numel(depth),1)=0;
for i=1:numel(depth)
    if i==1
        sigma_bulk_no_w(i,1)=(depth(i,1)*g*bulk(i,1));
        
    else
        
        sigma_bulk_no_w(i,1)= depth(i,1)*g*(mean(bulk(1:i,1))-water_dens);
    
    end
end
clear i

%convert to MPa
sigma_w_mpa=sigma_w./1e6;
hydrostatic_mpa=hydrostatic./1e6;

sigma_bulk_no_w_mpa=sigma_bulk_no_w./1e6;

% temperature calculations
temp_all=[34.6,21,0.15;
    63,40.7,2.88;
    120,89.5,4.96;
    148.5,111.5,5.95;
    189.1,149.3,7.7;
    234.1,182.8,10.28];
temp_depth=temp_all(:,1);
thermal_res=temp_all(:,2);
sampled_temp=temp_all(:,3);

temp_depth=temp_depth(~isnan(temp_depth));
sampled_temp=sampled_temp(~isnan(sampled_temp));

used_temp=sampled_temp(2:5);
used_temp_depth=temp_depth(2:5);

ft=fitlm(used_temp,used_temp_depth,'linear');

slope=table2array(ft.Coefficients(2,1));
intercept=table2array(ft.Coefficients(1,1));

%calculate temperature gradient
contact_depths = [110.5,222,510.96,848.45];
temp_estimate_depths=[0 509.82 550 600 650 700 750 800 848.45 1000 1050];
temp_estimate=(temp_estimate_depths-intercept)./slope;
temp_grad=temp_estimate(10)-temp_estimate(1);

% plot
plt1=figure('units','normalized','outerposition',[0 0 0.9 0.7]);
subplot(3,4,[1,9])
% h1=scatter(mad,depth,20);
% hold on 
h2=scatter(bulk,depth,20);
%pbaspect([1 3 1])
ylim([0 850])
set(gca,'YDir','reverse','XAxisLocation','top')
xlabel('Density (kg m^{-3})')
ylabel('Depth (mbsf)')
set(gca,'FontSize',fd_fontsize)
box off
hold on
for i = 1:numel(contact_depths)
    L = plot(xlim,[contact_depths(i),contact_depths(i)],'Color',[0.8,0.8,0.8],'LineWidth',2);
    uistack(L,'bottom')
end
hold off

subplot(3,4,[2,10])
h3=plot(sigma_bulk_no_w_mpa,depth,...
    hydrostatic_mpa,depth);
ylim([0 850])
xlim([0 45])
set(gca,'YDir','reverse','XAxisLocation','top')
xlabel('Stress (MPa)')
box off
lgd3=legend('Eff. vertical stress',...
    'Hydrostatic stress',...
    'Location','northwest');
lgd3.FontSize=fd_fontsize;
set(gca,'FontSize',fd_fontsize)
legend('boxoff')
hold on
legend('AutoUpdate','off')
for i = 1:numel(contact_depths)
    L = plot(xlim,[contact_depths(i),contact_depths(i)],'Color',[0.8,0.8,0.8],'LineWidth',2);
    uistack(L,'bottom')
end
hold off

subplot(3,4,[3,11])
plot(20,300,'Color','none')
hold on
h4=plot(sampled_temp,temp_depth,'o',used_temp,used_temp_depth,'o',temp_estimate,temp_estimate_depths,'-k');
ylim([0 850])
xlim([-1 50])
set(gca,'YDir','reverse','XAxisLocation','top')
set(h4,{'MarkerFaceColor'},get(h4,'Color'))
xlabel(['Temperature (' char(0176) 'C)'])
box off
lgd1=legend('Temperatures','Measurements',...
    'Used for gradient',...
    'Gradient',...
    'Location','northeast');
lgd1.FontSize=fd_fontsize;
legend('boxoff')
set(gca,'FontSize',fd_fontsize)
text(temp_estimate(6)+2,temp_estimate_depths(6),['R^2 = ' num2str(ft.Rsquared.Ordinary)],'fontsize',fd_fontsize)
text(0.2,820,['\deltaT/\deltaz = ' num2str(temp_grad,3) '\circC km^{-1}'],'fontsize',fd_fontsize)
legend('AutoUpdate','off')
for i = 1:numel(contact_depths)
    L = plot(xlim,[contact_depths(i),contact_depths(i)],'Color',[0.8,0.8,0.8],'LineWidth',2);
    uistack(L,'bottom')
end
hold off

%% Age Model 
%decide what to use for age model
%1=all, 2=plankton+paleomag, 3=plankton+nanofossil, 4=paleomag+nannofossil
%5=paleomag, 6=nannofossil, 7=plankton
option=1;

%averaging window in Ma
avg_win=0.25;

% import data
fossil_ages=xlsread('./data/U1520-D-T1 Datums.xlsx');
paleomag_age_depth=csvread('./data/paleomag_age_depth_no_repeats.csv');

%get age data
plankton_age=fossil_ages(:,16);
nano_age=fossil_ages(:,17);
paleomag_age=paleomag_age_depth(:,1);

%get fossil age errors
nano_age_pos_err=fossil_ages(:,13);
nano_age_neg_err=fossil_ages(:,12);
plankton_age_pos_err=fossil_ages(:,13);
plankton_age_neg_err=fossil_ages(:,12);

%get fossil depth errors
nano_depth_pos_err=fossil_ages(:,14);
nano_depth_neg_err=fossil_ages(:,15);
plankton_depth_pos_err=fossil_ages(:,14);
plankton_depth_neg_err=fossil_ages(:,15);


%sort errors into plankton and foram
nano_age_pos_err=nano_age_pos_err(isnan(nano_age)==0);
nano_age_neg_err=nano_age_neg_err(isnan(nano_age)==0);
plankton_age_pos_err=plankton_age_pos_err(isnan(plankton_age)==0);
plankton_age_neg_err=plankton_age_neg_err(isnan(plankton_age)==0);

nano_depth_pos_err=nano_depth_pos_err(isnan(nano_age)==0);
nano_depth_neg_err=nano_depth_neg_err(isnan(nano_age)==0);
plankton_depth_pos_err=plankton_depth_pos_err(isnan(plankton_age)==0);
plankton_depth_neg_err=plankton_depth_neg_err(isnan(plankton_age)==0);

%get depth data
plankton_depth=fossil_ages(:,11);
nano_depth=fossil_ages(:,11);
paleomag_depth=paleomag_age_depth(:,2);

%remove NaNs from fossil data
plankton_depth=plankton_depth(isnan(plankton_age)==0);
nano_depth=nano_depth(isnan(nano_age)==0);
plankton_age=plankton_age(isnan(plankton_age)==0);
nano_age=nano_age(isnan(nano_age)==0);


% set trigger for data

if option==1
    x=cat(1,plankton_age,nano_age,paleomag_age);
    y=cat(1,plankton_depth,nano_depth,paleomag_depth);
    age_str="Age model from all data";
elseif option==2
    x=cat(1,plankton_age,paleomag_age);
    y=cat(1,plankton_depth,paleomag_depth);  
        age_str="Age model from Forams and paleomag";
elseif option==3
    x=cat(1,plankton_age,nano_age);
    y=cat(1,plankton_depth,nano_depth);    
        age_str="Age model from Forams and nannofossils";
elseif option==4
    x=cat(1,nano_age,paleomag_age);
    y=cat(1,nano_depth,paleomag_depth);    
    age_str="Age model from nannofossils and paleomag";
elseif option==5
    x=paleomag_age;
    y=paleomag_depth;    
    age_str="Age model from paleomag";
elseif option==6
    x=nano_age;
    y=nano_depth;    
    age_str="Age model from nannofossils";
elseif option==7
    x=plankton_age;
    y=plankton_depth;   
    age_str="Age model from Forams";
else
    disp('Select a valid value for option.')
end


clear d1 m1 n1 c1 d1

x=sort(x);
y=sort(y);

% fitting age model
age = movmean(x,avg_win);

%linearise to 0
age(2:end+1)=age;
age(1)=0;
y(2:end+1)=y;
y(1)=0;

% plotting
% figure
subplot(3,4,[4,12])
errorbar(plankton_age,plankton_depth,...
    plankton_depth_neg_err,plankton_depth_pos_err,...
    plankton_age_neg_err,plankton_age_pos_err,'.','CapSize',4,'MarkerSize',10)
hold on
errorbar(nano_age,nano_depth,...
    nano_depth_neg_err,nano_depth_pos_err,...
    nano_age_neg_err,nano_age_pos_err,'.','CapSize',4,'MarkerSize',10)
plot(paleomag_age,paleomag_depth,'.','MarkerSize',10);
set(gca,'YDir','reverse','fontsize', fd_fontsize,'XAxisLocation','top')
box off
ylim([0 850])
xlim([-1.5 65])
xlabel('Sediment age (Ma)')
%ylabel('Depth (mbsf)')
hold on
plot(age,y,'-k','LineWidth',1.5)
lgd2=legend('Planktonic Foraminifera ages','Nannofossil ages',...
    'Palaeomag reversal ages',age_str,...
    'Location','northeast');
lgd2.FontSize= fd_fontsize-2;
legend('boxoff')
legend('AutoUpdate','off')
for i = 1:numel(contact_depths)
    L = plot(xlim,[contact_depths(i),contact_depths(i)],'Color',[0.8,0.8,0.8],'LineWidth',2);
    uistack(L,'bottom')
end

hold off

%% P-T model
%interpolate age of stress levels
%create unique points for age depth interpolation
y_short = y + rot90(linspace(0, 1, numel(y))*1E-10,3) ;
age_short = age + rot90(linspace(0, 1, numel(age))*1E-10,3) ;

depth_age=interp1(y_short,age_short,depth,'linear');
stress_age=interp1(y_short,age_short,depth,'linear');

sigma_bulk_no_w_mpa=flip(sigma_bulk_no_w_mpa);

% stress history
% require density with age and no NaNs
density_age = stress_age(~isnan(stress_age));
density_short=bulk;
density_short(numel(density_age)+1:end)=[];

% remove duplicates from both and make continously increasing 
% for interpolation
density_age=density_age + rot90(linspace(0, 1, numel(density_age))*1E-10,3) ;
density_short=density_short + rot90(linspace(0, 1, numel(density_short))*1E-10,3) ;

% remove duplicates from age and y
age_short = age + rot90(linspace(0, 1, numel(age))*1E-10,3) ;
y_short=y + rot90(linspace(0, 1, numel(y))*1E-10,3) ;

%age steps
stepsize=0.1;

% interpolate age of horizon
j=interp1(y_short,age_short,max(y_short),'linear');
k=interp1(y_short,age_short,848.4,'linear'); % bottom of Unit IV
l=interp1(y_short,age_short,509.82,'linear'); % top of Unit IV

% create time steps for modelling 
t_j=0.01:stepsize:j;
t_k=0.01:stepsize:k;
t_l=0.01:stepsize:l;

t_j(end+1)=j;
t_k(end+1)=k;
t_l(end+1)=l;

t_j=flip(t_j);
t_k=flip(t_k);
t_l=flip(t_l);

%append 0 for today
t_j(end+1)=0;
t_k(end+1)=0;
t_l(end+1)=0;

tightfig;

set(gcf,'Color','w')
%% mass loss and strain by lithology

%read major lithological horizons from simplified stratigraphy
major_litho=xlsread('./data/Simplified_domain_2_strat.xlsx');

lithos=[ "510.96-655.15 mbsf Calcareous mudstone"...
    "655.15-720.93 mbsf Marl" ...
    "720.93-738.68 mbsf Mudstone" ...
    "738.68-773.26 mbsf Marl" ...
    "773.26-781.12 mbsf Chalk" ...
    "781.12-797.99 mbsf Debris flow" ...
    "797.99-848.45 mbsf Chalk"];

top_of_unit=major_litho(1,1);
bottom_of_unit=major_litho(end,2);

% read core recovery per section
[core_recovery,core_rec_labels,core_rec_raw]=xlsread('./data/core_recovery_edited.xlsx');
top_cored=core_recovery(:,6);
bottom_cored=core_recovery(:,7);
recovery=core_recovery(:,8);
recovered_length=(bottom_cored-top_cored).*(recovery./100);

unit_boundary_depths(1:7,1)=major_litho(1:7,1);
unit_boundary_depths(8,1)=major_litho(7,2);
unit_boundary_depths(:,2)=unit_boundary_depths(:,1);

%calc mudstone cells are 1-19
%top marl cells are 20-27
%mudstone cells are 28-30
%bottom marl cells are 31-34
%top chalk cells are 35-36
%debris flow cells are 37-39
%bottom chalk cells are 40-45

% get recovered length per section
thicknesses(1)=sum(recovered_length(1:19));
thicknesses(2)=sum(recovered_length(20:27));
thicknesses(3)=sum(recovered_length(28:30));
thicknesses(4)=sum(recovered_length(31:34));
thicknesses(5)=sum(recovered_length(35:36));
thicknesses(6)=sum(recovered_length(37:39));
thicknesses(7)=sum(recovered_length(40:45));

thicknesses(8)=sum(recovered_length);
total_thickness=thicknesses(8);

thicknesses=rot90(thicknesses,3);

%get number of stylolites per lithology
styl_numbers=major_litho(:,5);
styl_numbers(end+1)=sum(styl_numbers);

%% load carb and plot that and phys props
carb_spread=xlsread('./data/Carbonates_12_2_2019_continuous.xlsx');
depth_carb=carb_spread(:,10);
orig_caco3=carb_spread(:,13);

avg_period=3;
carb=movmean(orig_caco3,avg_period);

%plot carb, porosity, clay content

%carb
figure('units','normalized','outerposition',[1 1 0.7 1]);

%porosity
porosity_mean=movmean(porosity,avg_period);
%clay
XRD=xlsread('/Users/harryleah/OneDrive - Cardiff University/IODP Expedition 375/U1520 Inputs/XRD/U1520-C-T2.MU_continuous.xlsx');
XRD_depths=XRD(:,5);
clay_content=XRD(:,21);
quartz_content=XRD(:,22);
fsp_content=XRD(:,23);

%caco3 zoom
subplot(2,3,[1,4])
plot(orig_caco3,depth_carb,'ob',...
    carb,depth_carb,'-k',...
    [0 100],unit_boundary_depths(1,:),':k',...
    [0 100],unit_boundary_depths(2,:),':k',...
    [0 100],unit_boundary_depths(3,:),':k',...
    [0 100],unit_boundary_depths(4,:),':k',...
    [0 100],unit_boundary_depths(5,:),':k',...
    [0 100],unit_boundary_depths(6,:),':k',...
    [0 100],unit_boundary_depths(7,:),':k',...
    [0 100],unit_boundary_depths(8,:),':k'...
    ,'MarkerSize',2,'LineWidth',1.2)
xlabel('CaCO_3 content (wt%)')
legend('Measurements',[num2str(avg_period) ' Period moving avg'],...
    'Location','best','Box','off')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out')
box off
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([10 100])

%porosity zoom
subplot(2,3,[2,5])
plot(porosity,depth,'ob',...
    porosity_mean,depth,'-k',...
    [0 100],unit_boundary_depths(1,:),':k',...
    [0 100],unit_boundary_depths(2,:),':k',...
    [0 100],unit_boundary_depths(3,:),':k',...
    [0 100],unit_boundary_depths(4,:),':k',...
    [0 100],unit_boundary_depths(5,:),':k',...
    [0 100],unit_boundary_depths(6,:),':k',...
    [0 100],unit_boundary_depths(7,:),':k',...
    [0 100],unit_boundary_depths(8,:),':k'...
    ,'MarkerSize',2,'LineWidth',1.2)
xlabel('Porosity (%)')
legend('Porosity measurements',[num2str(avg_period) ' Period moving avg'],...
    'Location','best','Box','off')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out','ytick',[],'ycolor',[1 1 1])
box off
ax = gca;                   
ax.YAxis.Visible = 'off'; 
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([10 70])

%clay zoom
subplot(2,3,[3,6])
plot(clay_content,XRD_depths,'-k',...
    quartz_content,XRD_depths,'-b',...
    fsp_content,XRD_depths,'-r',...
    [0 100],unit_boundary_depths(1,:),':k',...
    [0 100],unit_boundary_depths(2,:),':k',...
    [0 100],unit_boundary_depths(3,:),':k',...
    [0 100],unit_boundary_depths(4,:),':k',...
    [0 100],unit_boundary_depths(5,:),':k',...
    [0 100],unit_boundary_depths(6,:),':k',...
    [0 100],unit_boundary_depths(7,:),':k',...
    [0 100],unit_boundary_depths(8,:),':k'...
    ,'MarkerSize',5)
xlabel('Mineral content (wt%)')
legend('Clay content (wt%)',...
    'Quartz content (wt%)',...
    'Feldspar content (wt%)',...
    'Location','best','Box','off')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out','ytick',[],'ycolor',[1 1 1])
box off
ax = gca;                   
ax.YAxis.Visible = 'off'; 
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 60])

%% Calculate strain from stylolites and plot against cac03
%get cac03 values within units
mean_unit_carb_content(1)=mean(orig_caco3(depth_carb >= unit_boundary_depths(1,1) & depth_carb < unit_boundary_depths(2,1)));
mean_unit_carb_content(2)=mean(orig_caco3(depth_carb >= unit_boundary_depths(2,1) & depth_carb < unit_boundary_depths(3,1)));
mean_unit_carb_content(3)=mean(orig_caco3(depth_carb >= unit_boundary_depths(3,1) & depth_carb < unit_boundary_depths(4,1)));
mean_unit_carb_content(4)=mean(orig_caco3(depth_carb >= unit_boundary_depths(4,1) & depth_carb < unit_boundary_depths(5,1)));
mean_unit_carb_content(5)=mean(orig_caco3(depth_carb >= unit_boundary_depths(5,1) & depth_carb < unit_boundary_depths(6,1)));
mean_unit_carb_content(6)=mean(orig_caco3(depth_carb >= unit_boundary_depths(6,1) & depth_carb < unit_boundary_depths(7,1)));
mean_unit_carb_content(7)=mean(orig_caco3(depth_carb >= unit_boundary_depths(7,1) & depth_carb < unit_boundary_depths(8,1)));
mean_unit_carb_content(8)=mean(orig_caco3(depth_carb >= unit_boundary_depths(1,1) & depth_carb < unit_boundary_depths(8,1)));
%strain rate = strain/strain time 
run=1;
%assuming bulk strain on stylolite of 0.7
for strain=[0.1 0.7 0.9]
styl_width=50;%seam width in um
stretch=1-strain;%stretch within seam
seams_per_stylolite=7;

styl_width=styl_width/1e6;%convert stylolite width to m

styl_thickness=styl_width*seams_per_stylolite;%individual stylolite thickness in m

total_styl_thickness=(styl_width*seams_per_stylolite)*styl_numbers(8);% thickness of all styls in stack

styl_thicknesses=(styl_width*seams_per_stylolite).*styl_numbers;% thickness of styls in each litho

lo_styl=styl_thicknesses./stretch; %original thickness of what is now stylolites

undeformed_thicknesses=(thicknesses-styl_thicknesses)+lo_styl; %undeformed thickness of each lithology

whole_stack_stretches=(thicknesses)./(undeformed_thicknesses);
whole_stack_strains=1-whole_stack_stretches;
%get array of times in ma
strain_time=1:1:interp1(y_short,age_short,max(y_short),'linear');
%convert it to seconds
strain_time=strain_time.*(60*60*24*365.25*1e6);
%get strain rate
strain_rate=strain./strain_time;
%back to ma for plotting
strain_time=strain_time./(60*60*24*365.25*1e6);

%get unit and stack bases
base_unit_depths=major_litho(:,2);
base_unit_depths(8)=major_litho(end,2);
%get depth_ages
base_unit_ages=interp1(y_short,age_short,base_unit_depths,'linear','extrap');
%do for entire stack
strain_times=[
    linspace(0,base_unit_ages(1),100)
    linspace(0,base_unit_ages(2),100)
    linspace(0,base_unit_ages(3),100)
    linspace(0,base_unit_ages(4),100)
    linspace(0,base_unit_ages(5),100)
    linspace(0,base_unit_ages(6),100)
    linspace(0,base_unit_ages(7),100)
    linspace(0,base_unit_ages(8),100)
    ];
%convert to seconds
strain_times=strain_times.*(60*60*24*365.25*1e6);
%get strain rate
strain_rates=whole_stack_strains./strain_times;
%convert back to Ma for plotting
strain_times=strain_times./(60*60*24*365.25*1e6);

%strain at bottom of unit IV
strain_time_848=interp1(y_short,age_short,848.4,'linear');
%convert to seconds
strain_time_848=strain_time_848*(60*60*24*365.25*1e6);
%get rate
strain_rate_848=strain/strain_time_848;
%back to ma for plotting
strain_time_848=strain_time_848/(60*60*24*365.25*1e6);

%strain rate at top of Unit IV
strain_time_509=interp1(y_short,age_short,509.82,'linear');
%convert to seconds
strain_time_509=strain_time_509*(60*60*24*365.25*1e6);
%get rate
strain_rate_509=strain/strain_time_509;
%back to ma for plotting
strain_time_509=strain_time_509/(60*60*24*365.25*1e6);

litho_names=[ "Calcareous mudstone"...
    "Upper marl" ...
    "Mudstone" ...
    "Lower marl" ...
    "Upper chalk" ...
    "Debris flow" ...
    "Lower chalk"];

styl_depthz=load('styl_depths.mat');
styl_depths(1,:)=styl_depthz.styl_depth;


for i=1:size(recovery)
    mean_section_carbs(i)=nanmean(...
        orig_caco3(depth_carb >= top_cored(i) &...
        depth_carb <= bottom_cored(i)));
    section_length_recovered(i)=(bottom_cored(i)-top_cored(i))*(recovery(i)/100);
    section_styl_count(i)=numel(styl_depths(...
        styl_depths > top_cored(i) &...
        styl_depths < bottom_cored(i)));
    section_styl_thickness(i)=section_styl_count(i)*(styl_width*seams_per_stylolite);
    lo_section_styl(i)=section_styl_thickness(i)/stretch;
    undeformed_section_thicknesses(i)=...
        (section_length_recovered(i)-section_styl_thickness(i))...
        +lo_section_styl(i);
    section_stretches(i)=(section_length_recovered(i))/...
        (undeformed_section_thicknesses(i));
    if section_styl_count(i)==0
        section_strains(i)=0;
    else
        section_strains(i)=1-section_stretches(i);
    end
end

if run==1
    strain_01_litho=whole_stack_strains;
    strain_01_sections=section_strains;
elseif run==2
    strain_06_litho=whole_stack_strains;
    strain_06_sections=section_strains;
elseif run==3
    strain_09_litho=whole_stack_strains;
    strain_09_sections=section_strains;
end

run=run+1;

end

full_length_init=sum(undeformed_section_thicknesses);
full_length_final=sum(section_length_recovered);
full_length_change=full_length_init-full_length_final;
full_strain=full_length_change/full_length_init;

clear run

whole_stack_strains=strain_06_litho;
section_strains=strain_06_sections;

plotting_strain_ranges_litho=[strain_01_litho strain_09_litho];
plotting_strain_x_litho=[rot90(mean_unit_carb_content,3) rot90(mean_unit_carb_content,3)];

figure
plotting_strain_ranges_section=[rot90(strain_01_sections,3) rot90(strain_09_sections,3)];
plotting_strain_x_section=[rot90(mean_section_carbs,3) rot90(mean_section_carbs,3)];

cmp=flip(parula(length(plotting_strain_x_section)));
colormap(cmp)
%cstr2=subplot(4,5,[14,20]);
for i=1:length(plotting_strain_x_section)
plot(plotting_strain_x_section(i,:),plotting_strain_ranges_section(i,:),...
    '-','Color',cmp(i,:),'LineWidth',1.5)
hold on
end
colormap(gca,cmp)
scatter(mean_section_carbs,section_strains,25,1:1:numel(recovery),'filled',...
    'MarkerEdgeColor','k','LineWidth',0.3)
set(gca,'fontsize',15,'TickDir','out')
pbaspect([1 1 1])
%title('Strain vs CaCO_3 content')
ylim([0 0.025])
str_xlim=xlim;
str_ylim=ylim;
box off
c=colorbar('eastoutside');
col_ticks(1)=min(c.Ticks);
col_ticks(2)=max(c.Ticks);
c.Ticks=[1 numel(recovery)];
c.TickLabels={'Top of core','Bottom of core'};
c.Direction='reverse';

exponential_type=1;

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( mean_section_carbs, section_strains );
% Set up fittype and options.
if exponential_type==2
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Normalize = 'on';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.00182870714642907 1.9871194499651 -0.0012841823598624 1.40046723895794];
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
elseif exponential_type==1
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Normalize = 'on';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.0025889632817082 1.30697051725559];
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end
legend('off')
% Plot fit with data.
xdata=linspace(str_xlim(1),str_xlim(2),200);
ydata=fitresult(xdata);
plot(xdata,ydata,'-k','LineWidth',1.5)
xlim(str_xlim)
ylim(str_ylim)
xlabel('Mean CaCO_3 content within core section (wt%)')
ylabel('Stylolite-hosted strain within core section')
text(str_xlim(1)*1.1,str_ylim(2)*0.82,['R^2 = ' num2str(gof.rsquare)],'FontSize',15)
if exponential_type==2
    text(str_xlim(1)*1.1,str_ylim(2)*0.9,['y = ' num2str(fitresult.a,3) 'e^{' num2str(fitresult.b,3) 'x} + ' num2str(fitresult.c,3) 'e^{' num2str(fitresult.d,3) 'x}'],'FontSize',15)
elseif exponential_type==1
    text(str_xlim(1)*1.1,str_ylim(2)*0.9,['y = ' num2str(fitresult.a,3) 'e^{' num2str(fitresult.b,3) 'x}'],'FontSize',15)
end
box off
hold off

%% porosity analysis
clear porosity_1 porosity_2
% gaussian width for smoothing
gauss = 4;

%image rotation
rotation_1=1;
rotation_2=10;

%bounds for sections in image 1
upper_seam_bnd_1=400;
lower_seam_bound_1=515;
%seam_crack_bound_1=403;

%bounds for sections in image 2
upper_seam_bnd_2=270;
lower_seam_bound_2=520;
seam_crack_bound_2=448;

% load images
%22R3W image
I1=imread('./data/burrow_montage_8.tiff');
I1=imrotate(I1,rotation_1);
I1=imcrop(I1,[1648 38 (1648-824) 1075]);

figure('Units','normalized','OuterPosition',[0 0 0.7 0.8],'Color',[1 1 1])
colormap(gray)
ax1=axes('Position',[0.05 0.1 0.25 0.5],'FontSize',15);
im1=imagesc(imadjust(I1));
axis equal
xticks([])
yticks([])
set(gca,'XColor','none','YColor','none','Color','none','FontSize',15)
pos=ax1.Position;
title('20R4W')

I1 = imgaussfilt(I1,gauss);

%19R1W1 image
I2=imread('./data/19R1W01-scale.tif');
I2=imrotate(I2,-rotation_2);
I2=imcrop(I2,[152 171 (776-152) (1033-171)]);

ax3=axes('Position',[0.5 0.1 0.25 0.5],'FontSize',15);
imagesc(I2)
axis equal
xticks([])
yticks([])
set(gca,'XColor','none','YColor','none','Color','none','FontSize',15)
title('19R1W')

I2 = imgaussfilt(I2,gauss);

%pre-allocate for speed
counts_1(1:length(I1(:,1)))=NaN;
counts_2(1:length(I2))=NaN;
porosity_1(1:length(I1(:,1)))=NaN;
porosity_2(1:length(I2))=NaN;
error_1(1:length(I1))=NaN;
error_2(1:length(I2))=NaN;

%threshold for porosity 0-255
por_thresh_1=graythresh(I1)*255;
por_thresh_2=graythresh(I2)*255;

%get porosity by row
for i=1:length(I1(:,1))
    counts_1(i)=sum(I1(i,:)<=por_thresh_1);
    porosity_1(i)=counts_1(i)/numel(I1(i,:));
    error_1(i) = sqrt(((std(porosity_1(i))/mean(porosity_1(i)))^2)+1);
end

for i=1:length(I2(:,1))
    counts_2(i)=sum(I2(i,:)<=por_thresh_2);
    porosity_2(i)=counts_2(i)/numel(I2(i,:));
    error_2(i) = sqrt(((std(porosity_2(i))/mean(porosity_2(i)))^2)+1);
end

%get below, above, and seam means with and without crack
below_mean_1=mean(porosity_1(lower_seam_bound_1:end),'omitnan');
above_mean_1=mean(porosity_1(1:upper_seam_bnd_1),'omitnan');
seam_mean_w_crack_1=mean(porosity_1(upper_seam_bnd_1:lower_seam_bound_1),'omitnan');
%seam_mean_wout_crack_1=mean(porosity_1(seam_crack_bound_1:lower_seam_bound_1),'omitnan');

below_mean_2=mean(porosity_2(lower_seam_bound_2:end),'omitnan');
above_mean_2=mean(porosity_2(1:upper_seam_bnd_2),'omitnan');
seam_mean_w_crack_2=mean(porosity_2(upper_seam_bnd_2:lower_seam_bound_2),'omitnan');
seam_mean_wout_crack_2=mean(porosity_2(upper_seam_bnd_2:seam_crack_bound_2),'omitnan');

sample_19r1w_depth=805.51;
sample_20r4w_depth=818.43;
darkness_of_grey=0.9;

sample_19r1w_depth_range(1)=max(depth(depth<sample_19r1w_depth));
sample_19r1w_depth_range(2)=min(depth(depth>sample_19r1w_depth));

sample_20r4w_depth_range(1)=max(depth(depth<sample_20r4w_depth));
sample_20r4w_depth_range(2)=min(depth(depth>sample_20r4w_depth));

sample_19r1w_por_range(1)=porosity(depth==sample_19r1w_depth_range(1))/100;
sample_19r1w_por_range(2)=porosity(depth==sample_19r1w_depth_range(2))/100;

sample_20r4w_por_range(1)=(porosity(depth==sample_20r4w_depth_range(1))/100);
sample_20r4w_por_range(2)=(porosity(depth==sample_20r4w_depth_range(2))/100);

pore_plot_20r4w_x=[sample_20r4w_por_range(1) sample_20r4w_por_range(2) sample_20r4w_por_range(2) sample_20r4w_por_range(1)];
pore_plot_20r4w_y=[0 0 length(I1(:,1)) length(I1(:,1))];

pore_plot_19r1w_x=[sample_19r1w_por_range(1) sample_19r1w_por_range(2) sample_19r1w_por_range(2) sample_19r1w_por_range(1)];
pore_plot_19r1w_y=[0 0 length(I2(:,1)) length(I2(:,1))];


ax2=axes('Position',[0.3 0.1 0.15 0.5],'FontSize',15);
fill(pore_plot_20r4w_x,pore_plot_20r4w_y,[darkness_of_grey darkness_of_grey darkness_of_grey],'LineStyle','none')
hold on
plot(porosity_1,1:length(I1(:,1)),'-k',...
    [below_mean_1 below_mean_1],[lower_seam_bound_1 length(I1(:,1))],'-b',...
    [above_mean_1 above_mean_1],[1 upper_seam_bnd_1],'-b',...
    [seam_mean_w_crack_1 seam_mean_w_crack_1],[upper_seam_bnd_1 lower_seam_bound_1],'-b',...
    'LineWidth',1)
xlabel('Porosity by pixel row')
set(gca,'YColor','none','YDir','reverse','FontSize',15,'Color','none')
yticks([])
xlim([0 1])
xticks([0 0.25 0.5 0.75])
ylim([0 length(I1(:,1))])
title(['20R4W porosity (k<' num2str(por_thresh_1) ')'])
text(0.3,mean([upper_seam_bnd_1,lower_seam_bound_1]),...
    ['\Delta\Phi = ' num2str(mean([above_mean_1 below_mean_1])-seam_mean_w_crack_1,2) ],...
    'fontsize',15)
hold off

ax4=axes('Position',[0.8 0.1 0.15 0.5],'FontSize',15);
fill(pore_plot_19r1w_x,pore_plot_19r1w_y,[darkness_of_grey darkness_of_grey darkness_of_grey],'LineStyle','none')
hold on
plot(porosity_2,1:length(I2(:,1)),'-k',...
    [below_mean_2 below_mean_2],[lower_seam_bound_2 length(I2(:,1))],'-b',...
    [above_mean_2 above_mean_2],[1 upper_seam_bnd_2],'-b',...
    [seam_mean_w_crack_2 seam_mean_w_crack_2],[upper_seam_bnd_2 lower_seam_bound_2],'-b',...
    [seam_mean_wout_crack_2 seam_mean_wout_crack_2],[upper_seam_bnd_2 seam_crack_bound_2],'-b',...
    'LineWidth',1)
xlabel('Porosity by pixel row')
set(gca,'YColor','none','YDir','reverse','FontSize',15)
yticks([])
xlim([0 1])
xticks([0 0.25 0.5 0.75])
ylim([0 length(I2(:,1))])
title(['19R1W porosity (k<' num2str(por_thresh_2) ')'])
text(-0.2,0.87*mean([upper_seam_bnd_2,seam_crack_bound_2]),...
    ['\Delta\Phi = ' num2str(mean([above_mean_2 below_mean_2])-seam_mean_wout_crack_2,2) ],...
    'fontsize',15)
text(0.289,0.97*seam_crack_bound_2,...
    ['\Delta\Phi = ' num2str(mean([above_mean_2 below_mean_2])-seam_mean_w_crack_2,2) ],...
    'fontsize',15)
hold off

c1=colorbar('south');
c1.Position=[0.0655 0.07 0.22 0.02];
caxis([0 255])
c1.Ticks=[0 por_thresh_1 255];

c2=colorbar('south');
c2.Position=[0.502 0.07 0.245 0.02];
caxis([0 255])
c2.Ticks=[0 por_thresh_2 255];

%% Pressure solution model with time, then depth
% saving_size is the size at which to save figures
saving_size=1e-4;
% strain rate history using model of zhang and spiers
% loops for strain rate scaling with relative volume percents (use_rel_vol)
% then loops for grainsize (d, in m)
for use_rel_vol=[3 2]
for d=[1e-5 5e-5 1.5e-4 1e-4]
clear ps_sr_t_j t_j_depths T_t_j D_t_j cc_sol_t_j sigma_n_t_j stress_t_j hydrostat_t_j
clear t_j_depths_unit_iv_bounds T_t_j_unit_iv_bounds D_t_j_unit_iv_bounds D_t_j_unit_iv_bounds sigma_n_t_j_unit_iv_bounds stress_t_j_unit_iv_bounds ps_sr_t_j_unit_iv_bounds hydrostat_t_j_unit_iv_bounds

%using diffusion controlled model
phi(1,:)=porosity;%movmean(porosity,3);%initial porosity value
q=2*ceil(max(phi));%critical porosity
a_d=100;
S=1e-9;
omega_cc=3.693e-5;%molar volume of calcite
omega_qtz=(22.69./1e6);%molar volume of quartz
omega_fsp=(101.31/1e6);
omega_clay=(170.15/1e6);
R=8.314;%gas constant


%create slightly increasing values for itnerpolation
carb_unique = orig_caco3 + rot90(linspace(0, 1, numel(orig_caco3))*1E-10,3) ;
depth_carb_unique = depth_carb + rot90(linspace(0, 1, numel(depth_carb))*1E-10,3) ;
quartz_content_unique = quartz_content + rot90(linspace(0, 1, numel(quartz_content))*1E-10,3) ;
XRD_depths_unique = XRD_depths + rot90(linspace(0, 1, numel(XRD_depths))*1E-10,3) ;
fsp_content_unique = fsp_content + rot90(linspace(0, 1, numel(fsp_content))*1E-10,3) ;
clay_content_unique = clay_content + rot90(linspace(0, 1, numel(clay_content))*1E-10,3) ;

%normalise depths to standard depth measurement
carb_norm=interp1(depth_carb_unique,carb_unique,depth);
quartz_norm=interp1(XRD_depths_unique,quartz_content_unique,depth);
fsp_norm=interp1(XRD_depths_unique,fsp_content_unique,depth);
clay_norm=interp1(XRD_depths_unique,clay_content_unique,depth);

%calculate cc vol
cc_dens=(2.71*1e3);
cc_vol=((carb_norm./100).*cc_dens);
%calculate qtz vol
qtz_dens=(2.65.*1e3);
qtz_vol=((quartz_norm./100).*qtz_dens);
%calculate fsp vol
fsp_dens=3.01*1e3;
fsp_vol=((fsp_norm./100).*fsp_dens);
%calculate clay vol
sm_dens=2.35*1e3;
clay_vol=((clay_norm./100).*sm_dens);

%get relative vols
qtz_rel_vol=qtz_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
cc_rel_vol=cc_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
fsp_rel_vol=fsp_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
clay_rel_vol=clay_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);


%function to decide how to scale strain rate
% use_rel_vol=2;
% 0=no scaling
% 1=scaled solubilities by all minerals
% 2=scale strain rate by cc volume
% 3=Use soluble phase scaling of zubstov et al., 2004

%pre-allocate sizes
ps_sr_t_j(1:numel(t_j),1:numel(depth))=NaN;
t_j_depths(1:numel(t_j),1:numel(depth))=NaN;
T_t_j(1:numel(t_j),1:numel(depth))=NaN;
D_t_j(1:numel(t_j),1:numel(depth))=NaN;
cc_sol_t_j(1:numel(t_j),1:numel(depth))=NaN;
sigma_n_t_j(1:numel(t_j),1:numel(depth))=NaN;
stress_t_j(1:numel(t_j),1:numel(depth))=NaN;
ps_sr_t_j(1:numel(t_j),1:numel(depth))=NaN;
hydrostat_t_j(1:numel(t_j),1:numel(depth))=NaN;

t_j_depths_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
hydrostat_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
T_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
D_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
cc_sol_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
stress_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
sigma_n_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
ps_sr_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;     

w_stress(1:numel(t_j))=NaN;
w_stress_z_err(1:numel(t_j))=NaN;
w_stress_age_err(1:numel(t_j))=NaN;
t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
step_t=stepsize*(60*60*24*365.25*1e6);
L_t_j(1:numel(t_j),1:numel(depth))=NaN;
cum_t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
r(1:numel(t_j),1:numel(depth))=NaN;
x_q(1:numel(t_j),1:numel(depth))=NaN;
Ac(1:numel(t_j),1:numel(depth))=NaN;
fd(1:numel(t_j),1:numel(depth))=NaN;
L0=1;

%porosity change
use_por_change=0;
delta_phi=-0.2;
time_seam=interp1(y_short,age_short,805.51,'linear');%in Ma
por_seam_change=delta_phi/time_seam;

%for each time calculate a strain rate at each cell 
%onyl calculate if cell is at an age less than 
for t=1:numel(t_j)
    for z=1:numel(depth)
        if isnan(depth_age(z))==1
            ps_sr_t_j(t,z)=NaN;
            t_j_depths(t,z)=NaN;
            T_t_j(t,z)=NaN;
            D_t_j(t,z)=NaN;
            cc_sol_t_j(t,z)=NaN;
            sigma_n_t_j(t,z)=NaN;
            stress_t_j(t,z)=NaN;
        elseif t_j(t)>depth_age(z)
            ps_sr_t_j(t,z)=NaN;
            t_j_depths(t,z)=NaN;
            T_t_j(t,z)=NaN;
            D_t_j(t,z)=NaN;
            cc_sol_t_j(t,z)=NaN;
            sigma_n_t_j(t,z)=NaN;
            stress_t_j(t,z)=NaN;
        else
            t_j_depths(t,z)=depth(z)-interp1(age_short,y_short,t_j(t),'linear');
            T_t_j(t,z)=((t_j_depths(t,z)-intercept)./slope)+273.15;
            D_t_j(t,z)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j(t,z))));
            cc_sol_t_j(t,z)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j(t,z) + 2839.319./T_t_j(t,z) + 71.595.*log10(T_t_j(t,z)))).*omega_cc).*1000;
            sigma_n_t_j(t,z)= t_j_depths(t,z)*g*(mean(bulk(1:z))-water_dens);
            ps_sr_t_j(t,z)=(a_d.*((D_t_j(t,z).*S.*cc_sol_t_j(t,z))./(d^3)).*((sigma_n_t_j(t,z).*omega_cc)./(R.*T_t_j(t,z)))).*(1./((q-(2.*phi(z))).^2));
        end
    end
    %check ages for top of Unit IV
    if t_j(t)>interp1(depth,depth_age,unit_boundary_depths(1,1))
        t_j_depths_unit_iv_bounds(t,1)=NaN;
        hydrostat_t_j_unit_iv_bounds(t,1)=NaN;
        T_t_j_unit_iv_bounds(t,1)=NaN;
        D_t_j_unit_iv_bounds(t,1)=NaN;
        cc_sol_t_j_unit_iv_bounds(t,1)=NaN;
        stress_t_j_unit_iv_bounds(t,1)=NaN;
        sigma_n_t_j_unit_iv_bounds(t,1)=NaN;
        ps_sr_t_j_unit_iv_bounds(t,1)=NaN;
    else
        %calculate for top of unit IV
        t_j_depths_unit_iv_bounds(t,1)=unit_boundary_depths(1,1)-interp1(age_short,y_short,t_j(t),'linear');
        if t_j_depths_unit_iv_bounds(t,1)<0
            t_j_depths_unit_iv_bounds(t,1)=NaN;
        end
        hydrostat_t_j_unit_iv_bounds(t,1)=water_dens*g*t_j_depths_unit_iv_bounds(t,1);
        T_t_j_unit_iv_bounds(t,1)=((t_j_depths_unit_iv_bounds(t,1)-intercept)./slope)+273.15;
        D_t_j_unit_iv_bounds(t,1)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j_unit_iv_bounds(t,1))));
        cc_sol_t_j_unit_iv_bounds(t,1)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j_unit_iv_bounds(t,1) + 2839.319./T_t_j_unit_iv_bounds(t,1) + 71.595.*log10(T_t_j_unit_iv_bounds(t,1)))).*omega_cc).*1000;
        stress_t_j_unit_iv_bounds(t,1)=mean(bulk(depth > unit_boundary_depths(1,1)))*g*t_j_depths_unit_iv_bounds(t,1);
        sigma_n_t_j_unit_iv_bounds(t,1)=(mean(bulk(depth > unit_boundary_depths(1,1)))-water_dens)*g*t_j_depths_unit_iv_bounds(t,1);
        ps_sr_t_j_unit_iv_bounds(t,1)=(a_d.*((D_t_j_unit_iv_bounds(t,1).*S.*cc_sol_t_j_unit_iv_bounds(t,1))./(d^3)).*((sigma_n_t_j_unit_iv_bounds(t,1).*omega_cc)./(R.*T_t_j_unit_iv_bounds(t,1)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(1,1)))).^2));
    end
    %check ages for bottom of Unit IV
    if t_j(t)>interp1(depth,depth_age,unit_boundary_depths(end,1))
        t_j_depths_unit_iv_bounds(t,2)=NaN;
        hydrostat_t_j_unit_iv_bounds(t,2)=NaN;
        T_t_j_unit_iv_bounds(t,2)=NaN;
        D_t_j_unit_iv_bounds(t,2)=NaN;
        cc_sol_t_j_unit_iv_bounds(t,2)=NaN;
        stress_t_j_unit_iv_bounds(t,2)=NaN;
        sigma_n_t_j_unit_iv_bounds(t,2)=NaN;
        ps_sr_t_j_unit_iv_bounds(t,2)=NaN;
    else
        %calculate for bottom of Unit IV
        t_j_depths_unit_iv_bounds(t,2)=unit_boundary_depths(end,1)-interp1(age_short,y_short,t_j(t),'linear');
        if t_j_depths_unit_iv_bounds(t,2)<0
            t_j_depths_unit_iv_bounds(t,2)=NaN;
        end
        hydrostat_t_j_unit_iv_bounds(t,2)=water_dens*g*t_j_depths_unit_iv_bounds(t,2);
        T_t_j_unit_iv_bounds(t,2)=((t_j_depths_unit_iv_bounds(t,2)-intercept)./slope)+273.15;
        D_t_j_unit_iv_bounds(t,2)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j_unit_iv_bounds(t,2))));
        cc_sol_t_j_unit_iv_bounds(t,2)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j_unit_iv_bounds(t,2) + 2839.319./T_t_j_unit_iv_bounds(t,2) + 71.595.*log10(T_t_j_unit_iv_bounds(t,2)))).*omega_cc).*1000;
        stress_t_j_unit_iv_bounds(t,2)=mean(bulk(depth > unit_boundary_depths(end,1)))*g*t_j_depths_unit_iv_bounds(t,2);
        sigma_n_t_j_unit_iv_bounds(t,2)=(mean(bulk(depth > unit_boundary_depths(end,1)))-water_dens)*g*t_j_depths_unit_iv_bounds(t,2);
        ps_sr_t_j_unit_iv_bounds(t,2)=(a_d.*((D_t_j_unit_iv_bounds(t,2).*S.*cc_sol_t_j_unit_iv_bounds(t,2))./(d^3)).*((sigma_n_t_j_unit_iv_bounds(t,2).*omega_cc)./(R.*T_t_j_unit_iv_bounds(t,2)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(end,1)))).^2));
    end
end

if use_rel_vol==2
    for t=1:numel(t_j)
        for z=1:numel(depth)
            ps_sr_t_j(t,z)=ps_sr_t_j(t,z)*cc_rel_vol(z);
        end
        ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*interp1(depth,cc_rel_vol,unit_boundary_depths(1,1));
        ps_sr_t_j_unit_iv_bounds(t,2)=ps_sr_t_j_unit_iv_bounds(t,2)*interp1(depth,cc_rel_vol,unit_boundary_depths(end,1));
    end
elseif use_rel_vol==3
    for t=1:numel(t_j)
        for z=1:numel(depth)
            if cc_rel_vol<=0.45
                ps_sr_t_j(t,z)=ps_sr_t_j(t,z)*(cc_rel_vol(z)/0.45);
            elseif cc_rel_vol>=0.75
                ps_sr_t_j(t,z)=ps_sr_t_j(t,z)*((-0.8*cc_rel_vol(z))+1.6);
            end
        end
        cc_rel_vol_unit_iv_top = interp1(depth,cc_rel_vol,unit_boundary_depths(1,1));
        cc_rel_vol_unit_iv_bot = interp1(depth,cc_rel_vol,unit_boundary_depths(end,1));
        if cc_rel_vol_unit_iv_top<=0.45
            ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*(cc_rel_vol_unit_iv_top/0.45);
        elseif cc_rel_vol_unit_iv_top>=0.75
            ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*((-0.8*cc_rel_vol_unit_iv_top)+1.6);
        end
        if cc_rel_vol_unit_iv_bot<=0.45
            ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*(cc_rel_vol_unit_iv_bot/0.45);
        elseif cc_rel_vol_unit_iv_bot>=0.75
            ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*((-0.8*cc_rel_vol_unit_iv_bot)+1.6);
        end
    end
end

%% integrate with time to get strain
clear final_t_j_strain final_t_j_strain_unit_iv_bounds cum_t_j_strain cum_t_j_strain_unit_iv_bounds
% pre-allocate variables for speed
%first for full depth 
t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
step_t=stepsize*(60*60*24*365.25*1e6);
L_t_j(1:numel(t_j),1:numel(depth))=NaN;
cum_t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
%then for unit IV bounds
t_j_strain_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
L_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
cum_t_j_strain_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
L0=1;

for z=1:numel(depth)
    for t=1:length(t_j)
        t_j_strain(t,z)=-ps_sr_t_j(t,z)*step_t;
        if t==1
            L_t_j(t,z)=L0+(t_j_strain(t,z)*L0);
        elseif isnan(L_t_j(t-1,z))==1
            L_t_j(t,z)=L0+(t_j_strain(t,z)*L0);
        else
            L_t_j(t,z)=L_t_j(t-1,z)+(t_j_strain(t,z)*L_t_j(t-1,z));
        end
        cum_t_j_strain(t,z)=-(L_t_j(t,z)-L0)/L0;
    end
end

final_t_j_strain=cum_t_j_strain(end,:);

for z=1:2
    for t=1:length(t_j)
        t_j_strain_unit_iv_bounds(t,z)=-ps_sr_t_j_unit_iv_bounds(t,z)*step_t;
        if t==1
            L_unit_iv_bounds(t,z)=L0+(t_j_strain_unit_iv_bounds(t,z)*L0);
        elseif isnan(L_unit_iv_bounds(t-1,z))==1
            L_unit_iv_bounds(t,z)=L0+(t_j_strain_unit_iv_bounds(t,z)*L0);
        else
            L_unit_iv_bounds(t,z)=L_unit_iv_bounds(t-1,z)+(t_j_strain_unit_iv_bounds(t,z)*L_unit_iv_bounds(t-1,z));
        end
        cum_t_j_strain_unit_iv_bounds(t,z)=-(L_unit_iv_bounds(t,z)-L0)/L0;
    end
end


%% plot conditions and model plots

fd_fontsize = 17;

figure('Units','normalized','OuterPosition',[0.1 0 0.9 1],'Color','w');
subplot(3,2,1)
plot(t_j,t_j_depths_unit_iv_bounds(:,2),'-k',t_j,t_j_depths_unit_iv_bounds(:,1),'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
ylabel('Sediment depth (mbsf)')
legend('Bottom of Unit IV','Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize)
xlim([0 65])
ylim([0 900])
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'FontSize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 900])
set(gca,'box','off')

subplot(3,2,3)
plot(t_j,T_t_j_unit_iv_bounds(:,2)-273.15,'-k',t_j,T_t_j_unit_iv_bounds(:,1)-273.15,'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
ylabel('Temperature (\circC)')
legend('Bottom of Unit IV','Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize)
xlim([0 65])
ylim([0 35])
hold on
yyaxis right
plot(t_j,log10(cc_sol_t_j_unit_iv_bounds(:,2)),'-b','DisplayName','Calcite solubility, bottom of Unit IV','LineWidth',1)
ylabel('Log_{10} solubility (m^3/m^3)')
xlim([0 65])
hold off

subplot(3,2,5)
plot(t_j,stress_t_j_unit_iv_bounds(:,2)./1e6,'-k',...
    t_j,stress_t_j_unit_iv_bounds(:,1)./1e6,'-.k',...
    t_j,hydrostat_t_j_unit_iv_bounds(:,2)./1e6,'-r',...
    t_j,hydrostat_t_j_unit_iv_bounds(:,1)./1e6,'-.r','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
xlabel('Time before present (Myr)')
ylabel('Stress (MPa)')
set(gca,'box','off')
xlim([0 65])
lgd7=legend('Lithostatic, bottom of Unit IV',...
    'Lithostatic, top of Unit IV','Hydrostratic, bottom of Unit IV',...
    'Hydrostratic, top of Unit IV',...
    'Location','northwest','Box','off','fontsize',fd_fontsize-1);
lgd7.NumColumns = 2;
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 56])
hold off

subplot(3,2,2)
plot(t_j,sigma_n_t_j_unit_iv_bounds(:,2)./1e6,'-r',...
    t_j,sigma_n_t_j_unit_iv_bounds(:,1)./1e6,'-.r','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
ylabel({'Effective normal stress','(MPa)'})
xlim([0 65])
ylim([0 9])
leg01=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize);
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 9])
set(gca,'box','off')

subplot(3,2,4)
plot(t_j,log10(ps_sr_t_j_unit_iv_bounds(:,2)),'-k',t_j,log10(ps_sr_t_j_unit_iv_bounds(:,1)),'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
ylabel('Log_{10} strain rate (s^{-1})')
lgd08=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','northwest','Box','off','fontsize',fd_fontsize);
xlim([0 65])
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([-20 -16])
set(gca,'box','off')

subplot(3,2,6)
plot(t_j,cum_t_j_strain_unit_iv_bounds(:,2),'-k',t_j,cum_t_j_strain_unit_iv_bounds(:,1),'-.k','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
xlabel('Time before present (Myr)')
ylabel('Cumulative strain')
lgd09=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','northwest','Box','off',...
    'AutoUpdate','off','fontsize',fd_fontsize);
xlim([0 65])
hold on

%mini plot to add
mini_xmin = 0;
mini_xmax = 1.9;
mini_ymin = 0;
mini_ymax = 0.009;

plot([mini_xmax mini_xmin mini_xmin mini_xmax mini_xmax],...
    [mini_ymin mini_ymin mini_ymax mini_ymax mini_ymin],'r','LineWidth',1.5)
hold off

ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','fontsize',fd_fontsize,'YTickLabel',[])
ylim([0 0.015])
set(gca,'box','off')
mini_plot_ax = axes('Position',[0.835 0.125 0.065 0.065],...
    'Color','none');
plot(t_j,cum_t_j_strain_unit_iv_bounds(:,1),'-.k','LineWidth',1)
xlim([mini_xmin mini_xmax])
ylim([mini_ymin mini_ymax])
mini_plot_ax.YTickLabel=num2str(str2double(mini_plot_ax.YTickLabel).*(10^(-3)));
mini_plot_ax.XDir = 'reverse';
mini_plot_ax.FontSize =fd_fontsize;
mini_plot_ax.XAxisLocation = 'top';
mini_plot_ax.TickDir = 'out';
mini_plot_ax.XColor = 'red';
mini_plot_ax.YColor = 'red';
mini_plot_ax.LineWidth = 1.2;

tightfig;

if d==1.5e-4
    final_strain_150=final_t_j_strain;
end

% save as eps
if d==saving_size
    if use_rel_vol == 3
        export_fig '-depsc' 'conds_history_with_model_zubstov_2004' '-painters'
    else
        export_fig '-depsc' 'conds_history_with_model' '-painters'
    end
end


%% check with strains from stylolite
temp=final_t_j_strain;

%get mean strains for small grain size model
mean_unit_model_strain(1,1)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(2,1)));
mean_unit_model_strain(1,2)=mean(temp(depth >= unit_boundary_depths(2,1) & depth < unit_boundary_depths(3,1)));
mean_unit_model_strain(1,3)=mean(temp(depth >= unit_boundary_depths(3,1) & depth < unit_boundary_depths(4,1)));
mean_unit_model_strain(1,4)=mean(temp(depth >= unit_boundary_depths(4,1) & depth < unit_boundary_depths(5,1)));
mean_unit_model_strain(1,5)=mean(temp(depth >= unit_boundary_depths(5,1) & depth < unit_boundary_depths(6,1)));
mean_unit_model_strain(1,6)=mean(temp(depth >= unit_boundary_depths(6,1) & depth < unit_boundary_depths(7,1)));
mean_unit_model_strain(1,7)=mean(temp(depth >= unit_boundary_depths(7,1) & depth < unit_boundary_depths(8,1)));
mean_unit_model_strain(1,8)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(8,1)));

mean_unit_model_strain(2,:)=mean_unit_model_strain(1,:);

%for plotting out lithological boundaries
xlimss=[10^(-8) 10^(-3.75)];
if d==saving_size && use_rel_vol==2
figure('units','normalized','outerposition',[0 1 1 1]);
subplot(2,3,[1,4])
plot(final_t_j_strain,depth,'-k',...
    mean_unit_model_strain(1:2,8),[unit_boundary_depths(8,1) unit_boundary_depths(1,1)],':k',...
    xlimss,unit_boundary_depths(1,:),'--k',...
    xlimss,unit_boundary_depths(8,:),'--k')
xlabel('Strain from PS-mediated compaction model')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out')
ylim([min(depth) max(depth)])
legend(['d= ' num2str(d(1)*1e6) ' \mum'],'Mean strain for Unit IV')
legend('Box','off','Location','northeast')


subplot(2,3,[2,5])
plot(final_t_j_strain,depth,'-k',...
    mean_unit_model_strain(1:2,8),[unit_boundary_depths(8,1) unit_boundary_depths(1,1)],':k',...
    xlimss,unit_boundary_depths(1,:),'--k',...
    xlimss,unit_boundary_depths(2,:),'--k',...
    xlimss,unit_boundary_depths(3,:),'--k',...
    xlimss,unit_boundary_depths(4,:),'--k',...
    xlimss,unit_boundary_depths(5,:),'--k',...
    xlimss,unit_boundary_depths(6,:),'--k',...
    xlimss,unit_boundary_depths(7,:),'--k',...
    xlimss,unit_boundary_depths(8,:),'--k')
xlabel('Strain from PS-mediated compaction model')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out')
ylim([500 850])
legend(['d= ' num2str(d(1)*1e6) ' \mum'],'Mean strain for Unit IV')
legend('Box','off','Location','northeast')

subplot(2,3,[3,6])
plot(final_t_j_strain,depth,':k',...
    mean_unit_model_strain(1:2,1),[unit_boundary_depths(1,1) unit_boundary_depths(2,1)],'-k',...
    mean_unit_model_strain(1:2,2),[unit_boundary_depths(2,1) unit_boundary_depths(3,1)],'-k',...
    mean_unit_model_strain(1:2,3),[unit_boundary_depths(3,1) unit_boundary_depths(4,1)],'-k',...
    mean_unit_model_strain(1:2,4),[unit_boundary_depths(4,1) unit_boundary_depths(5,1)],'-k',...
    mean_unit_model_strain(1:2,5),[unit_boundary_depths(5,1) unit_boundary_depths(6,1)],'-k',...
    mean_unit_model_strain(1:2,6),[unit_boundary_depths(6,1) unit_boundary_depths(7,1)],'-k',...
    mean_unit_model_strain(1:2,7),[unit_boundary_depths(7,1) unit_boundary_depths(8,1)],'-k')
xlabel('Strain from PS-mediated compaction model')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'fontsize', 13,...
    'TickDir','out')
ylim([500 850])
    print '-depsc' 'model_final_strains' '-painters'
end
    final_for_later=final_t_j_strain;
    mean_unit_model_strain_later=mean_unit_model_strain;

%% compare strains to those from stylolites
%limits for axes
plot_points=[-4 -1];

mean_section_model_strains(1:numel(recovery))=NaN;
for i=1:numel(recovery)
    mean_section_model_strains(i)=mean(final_t_j_strain(depth >= top_cored(i) & depth <= bottom_cored(i)));
end

if d==saving_size && use_rel_vol==2
figure('units','normalized','outerposition',[1 0 0.8 0.6]);
subplot(1,2,1)
plot(log10(whole_stack_strains),log10(mean_unit_model_strain(1,:)),'.k','MarkerSize',25)
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','out')
legend(['d= ' num2str(d(1)*1e6) ' \mum'])
legend('Box','off','Location','southeast')
if use_rel_vol==1
    title('Comparing XRD-weighted modelled and stylolite based strains')
else
    title('Comparing modelled and stylolite strains')
end
hold on
label_shift=0.08;
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of calculated and modelled strains')
litho_stack_names=litho_names;
litho_stack_names(8)="Total Unit IV";
text(log10(whole_stack_strains)+label_shift,log10(mean_unit_model_strain(1,:)),litho_stack_names,'FontSize',15)%,'Rotation',45)
box off
hold off

plot_points=[-4.5 -1];

subplot(1,2,2)
scatter(log10(section_strains),log10(mean_section_model_strains),50,1:1:numel(recovery),'filled')
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','out')
xlim([-4.5 -1.5])
ylim([-4.5 -1.5])
colormap(flip(parula(45)))
c=colorbar('eastoutside');
col_ticks(1)=min(c.Ticks);
col_ticks(2)=max(c.Ticks);
c.Ticks=[1 numel(recovery)];
c.TickLabels={'Top of core','Bottom of core'};
c.Direction='reverse';
c.Label.FontSize=15;
hold on
plot(plot_points,plot_points,'-k')
hold off

    print '-depsc' 'calc_vs_model_strain_cc_scaled' '-painters'
end

%% trigger for grain sizes and plot
if d==1e-5
    litho_strains_10=mean_unit_model_strain(1,:);
    section_strains_10=mean_section_model_strains;
    close all
elseif d==5e-5
    litho_strains_50=mean_unit_model_strain(1,:);
    section_strains_50=mean_section_model_strains;
    close all
elseif d==1e-4
    litho_strains_100=mean_unit_model_strain(1,:);
    section_strains_100=mean_section_model_strains;
    close all
end


end
%% plot out comparison 

plot_points(1)=-4.5;
plot_points(2)=-0.5;

cmp=parula(10);
figure('units','normalized','outerposition',[1 0 0.8 0.6],'Color','w');
subplot(1,2,1)
plot(log10(whole_stack_strains),log10(litho_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(whole_stack_strains),log10(litho_strains_50),'.','Color',cmp(5,:),'MarkerSize',20)
plot(log10(whole_stack_strains),log10(litho_strains_100),'.','Color',cmp(8,:),'MarkerSize',20)
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','in')
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size')
legend('Box','off','Location','southeast','AutoUpdate','off')
hold on
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains')

litho_strains_plotting(:,1)=rot90(litho_strains_10,3);
litho_strains_plotting(:,2)=rot90(litho_strains_10,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(1,:))
end

litho_strains_plotting(:,1)=rot90(litho_strains_50,3);
litho_strains_plotting(:,2)=rot90(litho_strains_50,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(5,:))
end

litho_strains_plotting(:,1)=rot90(litho_strains_100,3);
litho_strains_plotting(:,2)=rot90(litho_strains_100,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(8,:))
end
hold off


subplot(1,2,2)
plot(log10(section_strains),log10(section_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(section_strains),log10(section_strains_50),'.','Color',cmp(5,:),'MarkerSize',20)
plot(log10(section_strains),log10(section_strains_100),'.','Color',cmp(8,:),'MarkerSize',20)
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','in')
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size')
legend('Box','off','Location','southeast','AutoUpdate','off')
hold on
label_shift=0.08;
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains')

section_strains_plotting(:,1)=rot90(section_strains_10,3);
section_strains_plotting(:,2)=rot90(section_strains_10,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(1,:))
end

section_strains_plotting(:,1)=rot90(section_strains_50,3);
section_strains_plotting(:,2)=rot90(section_strains_50,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(5,:))
end

section_strains_plotting(:,1)=rot90(section_strains_100,3);
section_strains_plotting(:,2)=rot90(section_strains_100,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(8,:))
end

hold off
tightfig;
if use_rel_vol == 3
    export_fig 'final_strain_comparison_var_grain_size_zubstov_2004.eps' '-painters' 
else
    export_fig 'final_strain_comparison_var_grain_size.eps' '-painters' 
end
%% plot images and final comparison
%create time-depth colour figures
cmp1=parula(20);
fd_plot=figure('Units','normalized','OuterPosition',[0.1 0 0.8 1],'Color','w');
axes('Position',[0.05 0.7 0.35 0.25])
p=pcolor(t_j,depth,transpose(sigma_n_t_j./1e6));
shading interp
colormap(gca,cmp1)
set(p,'AlphaData',transpose(~isnan(sigma_n_t_j)))
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
box off
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
c=colorbar('south');
c.Label.String = 'Effective Normal stress (MPa)';
c.Label.FontSize = fd_fontsize;
c.FontSize=fd_fontsize;
c.TickDirection='out';
c.Position=[0.075 0.92 0.25 0.025];
caxis([0 8])
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

cmp2=parula(20);

axes('Position',[0.455 0.7 0.35 0.25])
p=pcolor(t_j,depth,transpose(log10(cc_sol_t_j)));
shading interp
colormap(gca,cmp2)
set(p,'AlphaData',transpose(~isnan(sigma_n_t_j)))
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
box off
c=colorbar('south');
c.Label.String = 'Log_{10} Calcite solubility (m^3/m^3)';
c.Label.FontSize = fd_fontsize;
c.TickDirection='out';
c.Position=[0.47 0.92 0.25 0.025];
c.FontSize=fd_fontsize;
caxis([-5.715 -5.615])
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

cmp3=parula(14);

axes('Position',[0.05 0.06 0.35 0.575])
b=pcolor(t_j,depth,transpose(log10(ps_sr_t_j)));
shading interp
colormap(gca,cmp3)
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
box off
c=colorbar('south');
c.Label.String = 'Log_{10} modelled strain rate (s^{-1})';
c.Label.FontSize = fd_fontsize;
c.Position=[0.075 0.59 0.25 0.025];
c.TickDirection='out';
c.FontSize=fd_fontsize;
if d==5e-5
    caxis([-21 -14.5])
elseif d==1e-4
    caxis([-19 -15])
    colormap(gca,parula(16))
end
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

axes('Position',[0.455 0.06 0.35 0.575])
pre=pcolor(t_j,depth,transpose(cum_t_j_strain));
shading interp
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
if d==5e-5
    caxis([0 0.35])
colormap(gca,parula(14))
elseif d==1e-4
    caxis([0 0.15])
colormap(gca,parula(((length(0:0.025:0.15)*2)-2)*2))
end
box off
c=colorbar('south');
c.Label.String = 'Modelled cumulative strain';
c.Label.FontSize = fd_fontsize;
c.Position=[0.47 0.59 0.25 0.025];
c.TickDirection='out';
c.FontSize=fd_fontsize;
c.Ticks=0:0.025:0.15;
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

if d==saving_size
    if use_rel_vol == 3
        export_fig 'cum_model_strain_history_final_comp_zubstov_2004.png' -q101 -a1 '-r300'
    else
        export_fig 'cum_model_strain_history_final_comp.png' -q101 -a1 '-r300'
    end
end
% close

%% calculated concentrated clay per metre at given point from model
clear conc_dist conc_dist_removed clay_conc
conc_dist=1; %distance over which to calculate material removal in m
conc_dist_removed=(rot90(final_t_j_strain,3).*conc_dist);% distance of material removed in m (conc_dist./rot90((1-final_time_fut_strain),3))-conc_dist;
clay_conc=(clay_rel_vol.*conc_dist.*5.5.*5.5)...volume of clay in one metre of core (m3)
    ./((conc_dist-conc_dist_removed).*5.5.*5.5);%volume of remaining material per metre core (m3)
% diff_conc_orig=(clay_conc-clay_rel_vol);
% diff_conc_orig(diff_conc_orig<0)=0;
%get cumulative styl numbers per metre
max_depth=850;
min_depth=501;
cum_dist=1;

cum_styl(:,1:ceil(((max_depth-min_depth)/cum_dist)))=nan;

for i=1:ceil(((max_depth-min_depth)/cum_dist)+1)
    
    cum_styl(i)=sum(styl_depths<=(min_depth+(i*cum_dist))); 
    
end

cum_styl_grad(:,1:numel(cum_styl)-1)=nan;
cum_grad_depth(:,1:numel(cum_styl)-1)=nan;
depths=min_depth:1:max_depth;

for i=1:numel(cum_styl)
    if i<numel(cum_styl)
        cum_styl_grad(i)=cum_styl(i+1)-cum_styl(i);
    else
        cum_styl_grad(i)=0;
    end
    cum_grad_depth(i)=depths(i)+(cum_dist/2)+cum_dist;
    
    
end
cac03_p_m=interp1(depth_carb_unique,carb_unique,(depths-0.5));
strain=0.5;
styl_thickness_p_m=cum_styl_grad.*(styl_width*seams_per_stylolite);
stretch=1-strain;
lo_styl_p_m=styl_thickness_p_m./stretch;
undeformed_thicknesses_p_m=(conc_dist-styl_thickness_p_m)+lo_styl_p_m;
strain_p_m=1-(conc_dist./undeformed_thicknesses_p_m);
    
styl_dist_removed=((strain_p_m)*conc_dist);% distance of material removed in m (conc_dist./rot90((1-final_time_fut_strain),3))-conc_dist;
clay_styl=(interp1(depth,clay_rel_vol,(depths-0.5)).*conc_dist.*5.5.*5.5)...volume of clay in one metre of core (m3)
    ./((conc_dist-styl_dist_removed).*5.5.*5.5);

%% plot out comparison with final strains with depth

fd_fontsize = 17;

plot_points(1)=-5.5;
plot_points(2)=0;
cmp=parula(10);

figure('units','normalized','outerposition',[0 1 0.8 0.6],'Color','w');
subplot(2,5,[1 6])
plot(800,0.06,'Color','none')
hold on
plot(final_strain_150,depth,'-','Color',cmp(8,:),'LineWidth',1.5)
plot(final_t_j_strain,depth,'-','Color',cmp(6,:),'LineWidth',1.5)
plot(strain_p_m,(depths-0.5),'-','Color','k','LineWidth',1)
xlabel('Final strain per metre');

ylabel('Depth (mbsf)');

xticks([0 0.05  0.1 0.15])
set(gca,'YDir','reverse',...
    'XAxisLocation','bottom',...
    'fontsize', fd_fontsize,...
    'TickDir','in',...
    'YAxisLocation','left')
ylim([unit_boundary_depths(1,1) unit_boundary_depths(end,1)])
legend('Modelled',[num2str(150) ' \mum grain size'],[num2str(d*1e6) ' \mum grain size'],...
    'Observed','Box','Off','Location','northeast','FontSize',fd_fontsize)


subplot(2,5,[2 8])
plot(log10(whole_stack_strains),log10(litho_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(whole_stack_strains),log10(litho_strains_50),'.','Color',cmp(4,:),'MarkerSize',20)
plot(log10(whole_stack_strains),log10(litho_strains_100),'.','Color',cmp(6,:),'MarkerSize',20)
% Highlight total strain value for Unit IV
scatter(log10(whole_stack_strains(end)),log10(litho_strains_10(end)),40,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_50(end)),40,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_100(end)),40,'k','o')
xlim(plot_points)
ylim(plot_points)
yticks(-5:1:0)
xlabel('Log_{10} strain from stylolite mass loss and frequency');
T = ylabel('Log_{10} strain from pressure solution model');
T.VerticalAlignment = 'baseline';
T.Position(1) = T.Position(1)*1.035;
set(gca,'fontsize', fd_fontsize)
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size','Unit IV total')
legend('Box','off','Location','southeast','AutoUpdate','off','FontSize',fd_fontsize)
hold on
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains','LineWidth',1.5)

litho_strains_plotting(:,1)=rot90(litho_strains_10,3);
litho_strains_plotting(:,2)=rot90(litho_strains_10,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(1,:),'LineWidth',1)
end

litho_strains_plotting(:,1)=rot90(litho_strains_50,3);
litho_strains_plotting(:,2)=rot90(litho_strains_50,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(4,:),'LineWidth',1)
end

litho_strains_plotting(:,1)=rot90(litho_strains_100,3);
litho_strains_plotting(:,2)=rot90(litho_strains_100,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(6,:),'LineWidth',1)
end

line_lab = text(-3.5,-3.5,'Agreement of strains',...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',fd_fontsize);
line_lab.Rotation = 49;
hold off

subplot(2,5,[4 10])
plot(log10(section_strains),log10(section_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(section_strains),log10(section_strains_50),'.','Color',cmp(4,:),'MarkerSize',20)
plot(log10(section_strains),log10(section_strains_100),'.','Color',cmp(6,:),'MarkerSize',20)
xlim(plot_points)
ylim(plot_points)
xlabel('Log_{10} strain from stylolite mass loss and frequency');
set(gca,'fontsize', fd_fontsize)
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size')
legend('Box','off','Location','southeast','AutoUpdate','off','FontSize',fd_fontsize)
hold on
label_shift=0.08;
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains','LineWidth',1.5)

section_strains_plotting(:,1)=rot90(section_strains_10,3);
section_strains_plotting(:,2)=rot90(section_strains_10,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(1,:),'LineWidth',1)
end

section_strains_plotting(:,1)=rot90(section_strains_50,3);
section_strains_plotting(:,2)=rot90(section_strains_50,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(4,:),'LineWidth',1)
end

section_strains_plotting(:,1)=rot90(section_strains_100,3);
section_strains_plotting(:,2)=rot90(section_strains_100,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(6,:),'LineWidth',1)
end
yticks(-5:1:0)
line_lab = text(-3.5,-3.5,'Agreement of strains',...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',fd_fontsize);
line_lab.Rotation = 49;

hold off
tightfig;
if use_rel_vol ==3
    export_fig 'final_strain_comparison_with_depth_zubstov_2004.eps' '-painters' 
else
    export_fig 'final_strain_comparison_with_depth.eps' '-painters' 
end
%close
end
%% calculate into future
clear ps_sr_time_fut time_fut_depths T_time_fut cc_sol_time_fut sigma_n_time_fut stress_time_fut ps_sr_time_fut hydrostat_time_fut cum_time_fut_strain
final_t_j_strain_unit_iv_bounds=t_j_strain_unit_iv_bounds(end,:);
distance_step=100;  %in m
convergence_rate=40; %in mm/yr or km/myr
%import distances
thrust_dist=csvread('./data/distance to thrust.csv');
dist=-(thrust_dist(:,1)-max(thrust_dist(:,1)));
distance_to_toe=dist(2);
%use seismic constraints for distance to toe and sedimentation
seis_seafloor=csvread('./data/Seafloor_topUnitI.csv');
seis_seafloor_x=seis_seafloor(:,1);
seis_seafloor_y=seis_seafloor(:,2);
seis_seafloor_x=flip(seis_seafloor_x-min(seis_seafloor_x));

seis_top_U4=csvread('./data/Top Unit IV.csv');
seis_top_U4_x=seis_top_U4(:,1);
seis_top_U4_y=seis_top_U4(:,2);
seis_top_U4_x=flip(seis_top_U4_x-min(seis_top_U4_x));

distance_step=distance_step/1000;
modelled_distances=0:distance_step:distance_to_toe;
modelled_distances(end+1)=distance_to_toe; %in km

sedimentation=(interp1(seis_top_U4_x,seis_top_U4_y,modelled_distances)-interp1(seis_seafloor_x,seis_seafloor_y,modelled_distances));
sedimentation=sedimentation./(sedimentation(1)/unit_boundary_depths(1,1)); %in m
sea_level_change=1; %as a multiplier e.g. 1.1 causes a change of 10% each time step
%blank slate? 1=use zeros for current strain, 2=use final strain up to now
blank_slate=2;

%**** model*****

if blank_slate==2
    current_strain=final_t_j_strain;
    current_strain_unit_iv_bounds=final_t_j_strain_unit_iv_bounds;
else
    current_strain(1,1:numel(final_t_j_strain))=0;
    current_strain_unit_iv_bounds(size(final_t_j_strain_unit_iv_bounds))=0;
end
time_fut=(modelled_distances./convergence_rate); %in ma
sedimentation_rate=gradient(sedimentation)./gradient(time_fut);

%pre-allocate sizes
ps_sr_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
time_fut_depths(1:numel(time_fut),1:numel(depth))=NaN;
T_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
D_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
cc_sol_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
sigma_n_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
stress_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
ps_sr_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
hydrostat_time_fut(1:numel(time_fut),1:numel(depth))=NaN;

time_fut_depths_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
hydrostat_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
T_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
D_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
cc_sol_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
stress_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
sigma_n_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
ps_sr_time_fut_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;

sedimentation_diff=sedimentation-unit_boundary_depths(1,1);
sedimentation_diff(sedimentation_diff<0)=0;

for i=1:length(sedimentation_diff)
    if i>1
       if sedimentation_diff(i)<sedimentation_diff(i-1)
           sedimentation_diff(i)=sedimentation_diff(i-1);
       end
    end
end

w_stress(1:numel(time_fut))=NaN;
water_depth_fut=flip(interp1(seis_seafloor_x,seis_seafloor_y,modelled_distances)./(seis_seafloor_y(1)/water_depth));
%for each time calculate a strain rate at each cell 
for t=1:numel(time_fut)
    water_depth_fut(t)=water_depth_fut(t)*sea_level_change;
    for z=1:numel(depth)
            w_stress(t)=water_depth_fut(t)*water_dens*g;
            time_fut_depths(t,z)=depth(z)+(sedimentation_diff(t));
            hydrostatic=(time_fut_depths(t,z).*water_dens.*g)+w_stress(t);
            hydrostat_time_fut(t,z)=(water_depth_fut(t)+time_fut_depths(t,z))*g*water_dens;%
            T_time_fut(t,z)=((time_fut_depths(t,z)-intercept)./slope)+273.15;
            D_time_fut(t,z)=4.229e-08.*(exp(-15000./(8.3145.*T_time_fut(t,z))));
            cc_sol_time_fut(t,z)=(sqrt(10.^(-171.9065 - 0.077993.*T_time_fut(t,z) + 2839.319./T_time_fut(t,z) + 71.595.*log10(T_time_fut(t,z)))).*omega_cc).*1000;
            stress_time_fut(t,z)=(mean(bulk(1:z)*g*time_fut_depths(t,z)))+w_stress(t);
            sigma_n_time_fut(t,z)=stress_time_fut(t,z)-hydrostat_time_fut(t,z);
            ps_sr_time_fut(t,z)=(a_d.*((D_time_fut(t,z).*S.*cc_sol_time_fut(t,z))./(d^3)).*((sigma_n_time_fut(t,z).*omega_cc)./(R.*T_time_fut(t,z)))).*(1./((q-(2.*phi(z))).^2));

    end
    %calculate for top of unit IV
    time_fut_depths_unit_iv_bounds(t,1)=unit_boundary_depths(1,1)+(sedimentation_diff(t));
    hydrostat_time_fut_unit_iv_bounds(t,1)=(water_depth_fut(t)+time_fut_depths_unit_iv_bounds(t,1))*g*water_dens;%
    T_time_fut_unit_iv_bounds(t,1)=((time_fut_depths_unit_iv_bounds(t,1)-intercept)./slope)+273.15;
    D_time_fut_unit_iv_bounds(t,1)=4.229e-08.*(exp(-15000./(8.3145.*T_time_fut_unit_iv_bounds(t,1))));
    cc_sol_time_fut_unit_iv_bounds(t,1)=(sqrt(10.^(-171.9065 - 0.077993.*T_time_fut_unit_iv_bounds(t,1) + 2839.319./T_time_fut_unit_iv_bounds(t,1) + 71.595.*log10(T_time_fut_unit_iv_bounds(t,1)))).*omega_cc).*1000;
    stress_time_fut_unit_iv_bounds(t,1)=(mean(bulk(depth > unit_boundary_depths(1,1))*g*time_fut_depths_unit_iv_bounds(t,1)))+w_stress(t);
    sigma_n_time_fut_unit_iv_bounds(t,1)=stress_time_fut_unit_iv_bounds(t,1)-hydrostat_time_fut_unit_iv_bounds(t,1);
    ps_sr_time_fut_unit_iv_bounds(t,1)=(a_d.*((D_time_fut_unit_iv_bounds(t,1).*S.*cc_sol_time_fut_unit_iv_bounds(t,1))./(d^3)).*((sigma_n_time_fut_unit_iv_bounds(t,1).*omega_cc)./(R.*T_time_fut_unit_iv_bounds(t,1)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(1,1)))).^2));

    time_fut_depths_unit_iv_bounds(t,2)=unit_boundary_depths(end,1)+(sedimentation_diff(t));
    hydrostat_time_fut_unit_iv_bounds(t,2)=(water_depth_fut(t)+time_fut_depths_unit_iv_bounds(t,2))*g*water_dens;%
    T_time_fut_unit_iv_bounds(t,2)=((time_fut_depths_unit_iv_bounds(t,2)-intercept)./slope)+273.15;
    D_time_fut_unit_iv_bounds(t,2)=4.229e-08.*(exp(-15000./(8.3145.*T_time_fut_unit_iv_bounds(t,2))));
    cc_sol_time_fut_unit_iv_bounds(t,2)=(sqrt(10.^(-171.9065 - 0.077993.*T_time_fut_unit_iv_bounds(t,2) + 2839.319./T_time_fut_unit_iv_bounds(t,2) + 71.595.*log10(T_time_fut_unit_iv_bounds(t,2)))).*omega_cc).*1000;
    stress_time_fut_unit_iv_bounds(t,2)=(mean(bulk(depth > unit_boundary_depths(end,1))*g*time_fut_depths_unit_iv_bounds(t,2)))+w_stress(t);
    sigma_n_time_fut_unit_iv_bounds(t,2)=stress_time_fut_unit_iv_bounds(t,2)-hydrostat_time_fut_unit_iv_bounds(t,2);
    ps_sr_time_fut_unit_iv_bounds(t,2)=(a_d.*((D_time_fut_unit_iv_bounds(t,2).*S.*cc_sol_time_fut_unit_iv_bounds(t,2))./(d^3)).*((sigma_n_time_fut_unit_iv_bounds(t,2).*omega_cc)./(R.*T_time_fut_unit_iv_bounds(t,2)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(end,1)))).^2));
end

if use_rel_vol==2
    for t=1:numel(time_fut)
        for z=1:numel(depth)
            ps_sr_time_fut(t,z)=ps_sr_time_fut(t,z)*cc_rel_vol(z);
        end
        ps_sr_time_fut_unit_iv_bounds(t,1)=ps_sr_time_fut_unit_iv_bounds(t,1)*interp1(depth,cc_rel_vol,unit_boundary_depths(1,1));
        ps_sr_time_fut_unit_iv_bounds(t,2)=ps_sr_time_fut_unit_iv_bounds(t,2)*interp1(depth,cc_rel_vol,unit_boundary_depths(end,1));
    end
end

% calculate cumulative strains
step_fut=(distance_step/convergence_rate); %in Ma

% pre-allocate variables for speed
%first for full depth 
time_fut_strain(1:numel(time_fut),1:numel(depth))=NaN;
step_fut_s=step_fut*(60*60*24*365.25*1e6);
L_time_fut(1:numel(time_fut),1:numel(depth))=NaN;
cum_time_fut_strain(1:numel(time_fut),1:numel(depth))=NaN;
%then for unit IV bounds
time_fut_strain_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;
L_unit_iv_bounds_time_fut(1:numel(time_fut),1:2)=NaN;
cum_time_fut_strain_unit_iv_bounds(1:numel(time_fut),1:2)=NaN;

if blank_slate==2
    L=L_t_j(end,:);
else
    L(1:length(L_t_j(end,:)))=1;
end
L0=1;

for z=1:numel(depth)
    for t=1:length(time_fut)
        time_fut_strain(t,z)=-ps_sr_time_fut(t,z)*step_fut_s;
        if t==1
            L_time_fut(t,z)=L(z)+(time_fut_strain(t,z)*L(z));
        elseif isnan(L_time_fut(t-1,z))==1
            L_time_fut(t,z)=L(z)+(time_fut_strain(t,z)*L(z));
        else
            L_time_fut(t,z)=L_time_fut(t-1,z)+(time_fut_strain(t,z)*L_time_fut(t-1,z));
        end
        cum_time_fut_strain(t,z)=(-(L_time_fut(t,z)-L0)/L0);
    end
end

final_time_fut_strain=cum_time_fut_strain(end,:);

clear L0

L0=L_unit_iv_bounds(end,:);
for z=1:2
    for t=1:length(time_fut)
        time_fut_strain_unit_iv_bounds(t,z)=-ps_sr_time_fut_unit_iv_bounds(t,z)*step_fut_s;
        if t==1
            L_unit_iv_bounds_time_fut(t,z)=L0(z)+(time_fut_strain_unit_iv_bounds(t,z)*L0(z));
        elseif isnan(L_unit_iv_bounds_time_fut(t-1,z))==1
            L_unit_iv_bounds_time_fut(t,z)=L0(z)+(time_fut_strain_unit_iv_bounds(t,z)*L0(z));
        else
            L_unit_iv_bounds_time_fut(t,z)=L_unit_iv_bounds_time_fut(t-1,z)+(time_fut_strain_unit_iv_bounds(t,z)*L_unit_iv_bounds_time_fut(t-1,z));
        end
        cum_time_fut_strain_unit_iv_bounds(t,z)=(-(L_unit_iv_bounds_time_fut(t,z)-L0(z))/L0(z));
    end
end

temp=final_time_fut_strain;

%get mean strains for small grain size model
mean_unit_model_strain(1,1)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(2,1)));
mean_unit_model_strain(1,2)=mean(temp(depth >= unit_boundary_depths(2,1) & depth < unit_boundary_depths(3,1)));
mean_unit_model_strain(1,3)=mean(temp(depth >= unit_boundary_depths(3,1) & depth < unit_boundary_depths(4,1)));
mean_unit_model_strain(1,4)=mean(temp(depth >= unit_boundary_depths(4,1) & depth < unit_boundary_depths(5,1)));
mean_unit_model_strain(1,5)=mean(temp(depth >= unit_boundary_depths(5,1) & depth < unit_boundary_depths(6,1)));
mean_unit_model_strain(1,6)=mean(temp(depth >= unit_boundary_depths(6,1) & depth < unit_boundary_depths(7,1)));
mean_unit_model_strain(1,7)=mean(temp(depth >= unit_boundary_depths(7,1) & depth < unit_boundary_depths(8,1)));
mean_unit_model_strain(1,8)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(8,1)));

mean_unit_model_strain(2,:)=mean_unit_model_strain(1,:);

%% plot out future with distance
%images
fd_fontsize=18;
%create time-depth colour figures
figure('Units','normalized','OuterPosition',[0.1 0 0.8 1],'Color','w')
ax1=axes('Position',[0.05 0.69 0.4 0.25]);
p=pcolor(modelled_distances,depth,transpose(sigma_n_time_fut./1e6));
shading interp
set(p,'AlphaData',transpose(~isnan(sigma_n_time_fut)))
xlabel('Distance from Site U1520 (km)')
ylabel('Current depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','layer','top','fontsize',fd_fontsize)
box off
ylim([500 850])
colormap(gca,parula(18))
caxis([3 12])
c=colorbar('eastoutside');
c.Label.String = 'Effective Normal stress (MPa)';
c.Label.FontSize = fd_fontsize;
c.FontSize=fd_fontsize;
hold on
plot([0 max(modelled_distances)],unit_boundary_depths(2,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(3,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(4,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(5,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(6,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(7,:),'-k')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none','Nextplot','add','XLim',[0 max(time_fut)]);
yticks([])
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
xlabel('Time from present (Myr)')
hold off

ax1=axes('Position',[0.55 0.69 0.4 0.25]);
p=pcolor(modelled_distances,depth,transpose(log10(cc_sol_time_fut)));
shading interp
set(p,'AlphaData',transpose(~isnan(sigma_n_time_fut)))
xlabel('Distance from Site U1520 (km)')
ylabel('Current depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','layer','top','fontsize',fd_fontsize)
ylim([500 850])
box off
colormap(gca,parula(15))
caxis([-5.8 -5.65])
c=colorbar('eastoutside');
c.Label.String = 'Log_{10} Calcite solubility (m^3/m^3)';
c.Label.FontSize = fd_fontsize;
hold on
plot([0 max(modelled_distances)],unit_boundary_depths(2,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(3,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(4,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(5,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(6,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(7,:),'-k')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none','Nextplot','add','XLim',[0 max(time_fut)]);
yticks([])
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
xlabel('Time from present (Myr)')
hold off

ax1=axes('Position',[0.05 0.06 0.4 0.5]);
b=pcolor(modelled_distances,depth,transpose(log10(ps_sr_time_fut)));
shading interp
set(b,'AlphaData',transpose(~isnan(ps_sr_time_fut)))
xlabel('Distance from Site U1520 (km)')
ylabel('Current depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','layer','top','fontsize',fd_fontsize)
ylim([500 850])
box off
if d==5e-5
    colormap(gca,parula(20))
    caxis([-17 -15])
elseif d==1e-4
    caxis([-16 -14.75])
    colormap(gca,parula(10))
end
c=colorbar('eastoutside');
c.Label.String = 'Log_{10} modelled strain rate (s^{-1})';
c.Label.FontSize = fd_fontsize;
c.FontSize=fd_fontsize;
c.Ticks = -16:0.25:-14;
hold on
plot([0 max(modelled_distances)],unit_boundary_depths(2,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(3,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(4,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(5,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(6,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(7,:),'-k')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none','Nextplot','add','XLim',[0 max(time_fut)]);
yticks([])
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
xlabel('Time from present (Myr)')
hold off

ax1=axes('Position',[0.55 0.06 0.4 0.5]);
pre=pcolor(modelled_distances,depth,transpose(cum_time_fut_strain));
shading interp
set(pre,'AlphaData',transpose(~isnan(cum_time_fut_strain)))
xlabel('Distance from Site U1520 (km)')
ylabel('Current depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','layer','top','fontsize',fd_fontsize)
ylim([500 850])
box off
if d==5e-5
colormap(gca,parula(7))
caxis([0 0.35])
elseif d==1e-4
    caxis([0 0.175])
    colormap(gca,parula(14))
end
c=colorbar('eastoutside');
c.Label.String = 'Modelled cumulative strain';
c.Label.FontSize = fd_fontsize;
c.FontSize=fd_fontsize;
c.Ticks = 0:0.025:0.175;
hold on
plot([0 max(modelled_distances)],unit_boundary_depths(2,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(3,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(4,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(5,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(6,:),'-k',...
    [0 max(modelled_distances)],unit_boundary_depths(7,:),'-k')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Color','none','Nextplot','add','XLim',[0 max(time_fut)]);
yticks([])
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
xlabel('Time from present (Myr)')
hold off

if blank_slate==1
    export_fig 'cum_future_strain_dist_blank_slate.png' -q101 -a1'-r300'
elseif blank_slate==2
    export_fig 'cum_future_strain_dist.png' -q101 -a1 '-r300'
end

Time_until_sed_reaches_toe=[num2str(max(time_fut),3) ' Myr'];

%% plot calcite solubility with distance from U1520 down interface
%limit of depth
depth_lim=9;

temp_T=csvread('hik_temp_antriasian_2018.csv');
%define vars
dist_t_init=temp_T(:,1);   %distance from trench
plate_temp_init=temp_T(:,2);   %temperature

clear seis_seafloor
seis_seafloor=csvread('./data/Seafloor_topUnitI.csv');
seis_seafloor_x=seis_seafloor(:,1);
seis_seafloor_y=seis_seafloor(:,2);
seis_seafloor_x=flip(seis_seafloor_x-min(seis_seafloor_x));

geometry=csvread('./data/hik_depth_aintrasian.csv');
dist_d=geometry(:,1);
depth_init=geometry(:,2);

depth=depth_init(depth_init<depth_lim);
dist_d=dist_d(depth_init<depth_lim);

dist_t=dist_t_init(dist_t_init<max(dist_d));
plate_temp=plate_temp_init(dist_t_init<max(dist_d));
plate_temp=plate_temp+273.15;

seis_bot_U4=csvread('./data/Bottom Unit IV.csv');
seis_bot_U4_x=seis_bot_U4(:,1);
seis_bot_U4_y=seis_bot_U4(:,2);
seis_bot_U4_x=flip(seis_bot_U4_x-min(seis_bot_U4_x));
bot_u4_dist=0:0.1:max(seis_bot_U4_x);
bot_u4_depth=interp1(seis_bot_U4_x,seis_bot_U4_y,bot_u4_dist)-interp1(seis_seafloor_x,seis_seafloor_y,bot_u4_dist);


ft=fitlm(used_temp,used_temp_depth,'linear');
slope=table2array(ft.Coefficients(2,1));
intercept=table2array(ft.Coefficients(1,1));

T_bot_U4=(((bot_u4_depth.*1000)-intercept)./slope)+273.15;

thrust_dist=csvread('./data/distance to thrust.csv');
dist=-(thrust_dist(:,1)-max(thrust_dist(:,1)));
distance_to_toe=dist(2);
distance_to_intersect=dist(1);
decoll=csvread('./data/plate_interface.csv');
decoll_x=-(decoll(:,1)-max(decoll(:,1)));
decoll_t=interp1(dist_t,plate_temp,decoll_x);

plate_interface_sol=(sqrt(10.^(-171.9065 - 0.077993.*decoll_t + 2839.319./decoll_t + 71.595.*log10(decoll_t))).*omega_cc).*1000;
cc_sol_bot_U4=(sqrt(10.^(-171.9065 - 0.077993.*T_bot_U4 + 2839.319./T_bot_U4 + 71.595.*log10(T_bot_U4))).*omega_cc).*1000;

cc_sol_distances=cat(1,(decoll_x+distance_to_intersect),flip(rot90(bot_u4_dist,3)));
cc_sol_plate_interface=cat(1,plate_interface_sol,flip(rot90(cc_sol_bot_U4,3)));

cc_sol_distances(isnan(cc_sol_plate_interface)==1)=[];
cc_sol_plate_interface(isnan(cc_sol_plate_interface)==1)=[];

plate_sol=(sqrt(10.^(-171.9065 - 0.077993.*plate_temp + 2839.319./plate_temp + 71.595.*log10(plate_temp))).*omega_cc).*1000;

ftsz=18;
figure('Units','normalized','OuterPosition',[0 0 0.6 0.4])
plot(dist_t+distance_to_toe,log10(plate_sol),'-k',[20 30],[log10(cc_sol_t_j_unit_iv_bounds(end,2)) log10(cc_sol_t_j_unit_iv_bounds(end,2))],'--k')%,cc_sol_distances,log10(cc_sol_plate_interface),'--k')
ylabel('Log_{10} calcite solubility (m^3/m^3)')
xlabel('Distance from site U1520')
set(gca,'XDir','reverse','FontSize',18)
xlim([distance_to_toe 70])
hold on
text(30.5,log10(cc_sol_t_j_unit_iv_bounds(end,2)),'Bottom of Unit IV','FontSize',18,'HorizontalAlignment','right')
yyaxis right
plot(dist_t+distance_to_toe,plate_temp-273.15,'-r')
ylabel('Temperature (?C)')
set(gca,'YColor','r')
hold off

print '-depsc' 'sol_on_sed_interface.eps' '-painters'

%% get plot for calcite solubility and strain rate from model

strain_rate_at_interface_shear =(a_d.*(((4.229e-08.*(exp(-15000./(8.3145.*plate_temp)))).*S.*plate_sol)./(d^3)).*((2e6.*omega_cc)./(R.*plate_temp))).*(1./((q-(2.*35)).^2));
strain_rate_at_interface_actual =(a_d.*(((4.229e-08.*(exp(-15000./(8.3145.*plate_temp)))).*S.*plate_sol)./(d^3)).*((10e6.*omega_cc)./(R.*plate_temp))).*(1./((q-(2.*35)).^2));

ftsz=18;
figure('Units','normalized','OuterPosition',[0 0 0.6 0.4])
plot(dist_t+distance_to_toe,log10(plate_sol),'-k')
ylabel('Log_{10} calcite solubility (m^3/m^3)')
xlabel('Distance from site U1520')
set(gca,'XDir','reverse','FontSize',18)
xlim([distance_to_toe 70])
hold on
yyaxis right
plot(dist_t+distance_to_toe,log10(strain_rate_at_interface_shear),'--r',...
    dist_t+distance_to_toe,log10(strain_rate_at_interface_actual),'-r')
ylabel('Log_{10} strain rate (s^{-1})')
set(gca,'YColor','r')
legend('Calcite Solubility','Shear strain rate from PS at 2 MPa shear stress','Uniaxial strain rate from PS at 10 MPa normal stress','Location','northwest','Box','off')
hold off

%% Pressure solution model with time, then depth to determine possible overpressure

%strain rate history using model of zhang and spiers
close all
for d=[1e-5 5e-5 1.5e-4 1e-4]
%using diffusion controlled model

%decide if doing pore fluid factor analysis
use_pff = 1;
% 1 is true
% 0 is false

%pick a pore fluid factor
pff = 0.95;%1025/1900=0.5395
%agreement at 0.99 for 50 micron grain size
%agreement at 0.999 for 10 micron grain size

phi(1,:)=porosity;%initial porosity value
q=2*ceil(max(phi));%critical porosity
a_d=100;
S=1e-9;
% sigma_n=sigma_bulk_no_w-hydrostatic;%effective normal stress perpendicular to stylolite
omega_cc=3.693e-5;%molar volume of calcite
omega_qtz=(22.69./1e6);%molar volume of quartz
omega_fsp=(101.31/1e6);
omega_clay=(170.15/1e6);
R=8.31;%gas constant
depth=S1(:,9);      %depth in m

%create slightly increasing values for itnerpolation
carb_unique = orig_caco3 + rot90(linspace(0, 1, numel(orig_caco3))*1E-10,3) ;
depth_carb_unique = depth_carb + rot90(linspace(0, 1, numel(depth_carb))*1E-10,3) ;
quartz_content_unique = quartz_content + rot90(linspace(0, 1, numel(quartz_content))*1E-10,3) ;
XRD_depths_unique = XRD_depths + rot90(linspace(0, 1, numel(XRD_depths))*1E-10,3) ;
fsp_content_unique = fsp_content + rot90(linspace(0, 1, numel(fsp_content))*1E-10,3) ;
clay_content_unique = clay_content + rot90(linspace(0, 1, numel(clay_content))*1E-10,3) ;

%normalise depths to standard depth measurement
carb_norm=interp1(depth_carb_unique,carb_unique,depth);
quartz_norm=interp1(XRD_depths_unique,quartz_content_unique,depth);
fsp_norm=interp1(XRD_depths_unique,fsp_content_unique,depth);
clay_norm=interp1(XRD_depths_unique,clay_content_unique,depth);

%calculate cc vol
cc_dens=(2.71*1e3);
cc_vol=((carb_norm./100).*cc_dens);
%calculate qtz vol
qtz_dens=(2.65.*1e3);
qtz_vol=((quartz_norm./100).*qtz_dens);
%calculate fsp vol
fsp_dens=3.01*1e3;
fsp_vol=((fsp_norm./100).*fsp_dens);
%calculate clay vol
sm_dens=2.35*1e3;
clay_vol=((clay_norm./100).*sm_dens);

%get relative vols
qtz_rel_vol=qtz_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
cc_rel_vol=cc_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
fsp_rel_vol=fsp_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);
clay_rel_vol=clay_vol./(qtz_vol+cc_vol+fsp_vol+clay_vol);


%function to decide how to scale strain rate
use_rel_vol=2;
% 0=no scaling
% 1=scaled solubilities by all minerals
% 2=scale strain rate by cc volume

clear ps_sr_t_j t_j_depths T_t_j D_t_j cc_sol_t_j sigma_n_t_j stress_t_j hydrostat_t_j
%pre-allocate sizes
ps_sr_t_j(1:numel(t_j),1:numel(depth))=NaN;
t_j_depths(1:numel(t_j),1:numel(depth))=NaN;
T_t_j(1:numel(t_j),1:numel(depth))=NaN;
D_t_j(1:numel(t_j),1:numel(depth))=NaN;
cc_sol_t_j(1:numel(t_j),1:numel(depth))=NaN;
sigma_n_t_j(1:numel(t_j),1:numel(depth))=NaN;
stress_t_j(1:numel(t_j),1:numel(depth))=NaN;
hydrostat_t_j(1:numel(t_j),1:numel(depth))=NaN;

clear t_j_depths_unit_iv_bounds T_t_j_unit_iv_bounds D_t_j_unit_iv_bounds D_t_j_unit_iv_bounds sigma_n_t_j_unit_iv_bounds stress_t_j_unit_iv_bounds ps_sr_t_j_unit_iv_bounds hydrostat_t_j_unit_iv_bounds

t_j_depths_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
hydrostat_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
T_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
D_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
cc_sol_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
stress_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
sigma_n_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
ps_sr_t_j_unit_iv_bounds(1:numel(t_j),1:2)=NaN;     

w_stress(1:numel(t_j))=NaN;%get stress due to sea level
w_stress_z_err(1:numel(t_j))=NaN;
w_stress_age_err(1:numel(t_j))=NaN;
t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
step_t=stepsize*(60*60*24*365.25*1e6);
L_t_j(1:numel(t_j),1:numel(depth))=NaN;
cum_t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
r(1:numel(t_j),1:numel(depth))=NaN;
x_q(1:numel(t_j),1:numel(depth))=NaN;
Ac(1:numel(t_j),1:numel(depth))=NaN;
fd(1:numel(t_j),1:numel(depth))=NaN;
L0=1;

%porosity change
use_por_change=0;
delta_phi=-0.2;
time_seam=interp1(y_short,age_short,805.51,'linear');%in Ma
por_seam_change=delta_phi/time_seam;

%for each time calculate a strain rate at each cell 
%only calculate if cell is at an age less than 
for t=1:numel(t_j)
    for z=1:numel(depth)
        if isnan(depth_age(z))==1
            ps_sr_t_j(t,z)=NaN;
            t_j_depths(t,z)=NaN;
            T_t_j(t,z)=NaN;
            D_t_j(t,z)=NaN;
            cc_sol_t_j(t,z)=NaN;
            sigma_n_t_j(t,z)=NaN;
            stress_t_j(t,z)=NaN;
        elseif t_j(t)>depth_age(z)
            ps_sr_t_j(t,z)=NaN;
            t_j_depths(t,z)=NaN;
            T_t_j(t,z)=NaN;
            D_t_j(t,z)=NaN;
            cc_sol_t_j(t,z)=NaN;
            sigma_n_t_j(t,z)=NaN;
            stress_t_j(t,z)=NaN;
        else
            t_j_depths(t,z)=depth(z)-interp1(age_short,y_short,t_j(t),'linear');
            T_t_j(t,z)=((t_j_depths(t,z)-intercept)./slope)+273.15;
            D_t_j(t,z)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j(t,z))));
            cc_sol_t_j(t,z)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j(t,z) + 2839.319./T_t_j(t,z) + 71.595.*log10(T_t_j(t,z)))).*omega_cc).*1000;
            if use_pff == 1
                sigma_n_t_j(t,z)= t_j_depths(t,z)*g*(mean(bulk(1:z))-(mean(bulk(1:z))*pff));
            else
                sigma_n_t_j(t,z)= t_j_depths(t,z)*g*(mean(bulk(1:z))-water_dens);
            end
            ps_sr_t_j(t,z)=(a_d.*((D_t_j(t,z).*S.*cc_sol_t_j(t,z))./(d^3)).*((sigma_n_t_j(t,z).*omega_cc)./(R.*T_t_j(t,z)))).*(1./((q-(2.*phi(z))).^2));
        end
    end
    %check ages for top of Unit IV
    if t_j(t)>interp1(depth,depth_age,unit_boundary_depths(1,1))
        t_j_depths_unit_iv_bounds(t,1)=NaN;
        hydrostat_t_j_unit_iv_bounds(t,1)=NaN;
        T_t_j_unit_iv_bounds(t,1)=NaN;
        D_t_j_unit_iv_bounds(t,1)=NaN;
        cc_sol_t_j_unit_iv_bounds(t,1)=NaN;
        stress_t_j_unit_iv_bounds(t,1)=NaN;
        sigma_n_t_j_unit_iv_bounds(t,1)=NaN;
        ps_sr_t_j_unit_iv_bounds(t,1)=NaN;
    else
        %calculate for top of unit IV
        t_j_depths_unit_iv_bounds(t,1)=unit_boundary_depths(1,1)-interp1(age_short,y_short,t_j(t),'linear');
        if t_j_depths_unit_iv_bounds(t,1)<0
            t_j_depths_unit_iv_bounds(t,1)=NaN;
        end
        hydrostat_t_j_unit_iv_bounds(t,1)=water_dens*g*t_j_depths_unit_iv_bounds(t,1);
        T_t_j_unit_iv_bounds(t,1)=((t_j_depths_unit_iv_bounds(t,1)-intercept)./slope)+273.15;
        D_t_j_unit_iv_bounds(t,1)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j_unit_iv_bounds(t,1))));
        cc_sol_t_j_unit_iv_bounds(t,1)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j_unit_iv_bounds(t,1) + 2839.319./T_t_j_unit_iv_bounds(t,1) + 71.595.*log10(T_t_j_unit_iv_bounds(t,1)))).*omega_cc).*1000;
        stress_t_j_unit_iv_bounds(t,1)=mean(bulk(depth > unit_boundary_depths(1,1)))*g*t_j_depths_unit_iv_bounds(t,1);
        if use_pff == 1
            sigma_n_t_j_unit_iv_bounds(t,1)=(mean(bulk(depth > unit_boundary_depths(1,1)))-(mean(bulk(depth > unit_boundary_depths(1,1)))*pff))*g*t_j_depths_unit_iv_bounds(t,1);
        else
            sigma_n_t_j_unit_iv_bounds(t,1)=(mean(bulk(depth > unit_boundary_depths(1,1)))-water_dens)*g*t_j_depths_unit_iv_bounds(t,1);
        end
        ps_sr_t_j_unit_iv_bounds(t,1)=(a_d.*((D_t_j_unit_iv_bounds(t,1).*S.*cc_sol_t_j_unit_iv_bounds(t,1))./(d^3)).*((sigma_n_t_j_unit_iv_bounds(t,1).*omega_cc)./(R.*T_t_j_unit_iv_bounds(t,1)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(1,1)))).^2));
    end
    %check ages for bottom of Unit IV
    if t_j(t)>interp1(depth,depth_age,unit_boundary_depths(end,1))
        t_j_depths_unit_iv_bounds(t,2)=NaN;
        hydrostat_t_j_unit_iv_bounds(t,2)=NaN;
        T_t_j_unit_iv_bounds(t,2)=NaN;
        D_t_j_unit_iv_bounds(t,2)=NaN;
        cc_sol_t_j_unit_iv_bounds(t,2)=NaN;
        stress_t_j_unit_iv_bounds(t,2)=NaN;
        sigma_n_t_j_unit_iv_bounds(t,2)=NaN;
        ps_sr_t_j_unit_iv_bounds(t,2)=NaN;
    else
        %calculate for bottom of Unit IV
        t_j_depths_unit_iv_bounds(t,2)=unit_boundary_depths(end,1)-interp1(age_short,y_short,t_j(t),'linear');
        if t_j_depths_unit_iv_bounds(t,2)<0
            t_j_depths_unit_iv_bounds(t,2)=NaN;
        end
        hydrostat_t_j_unit_iv_bounds(t,2)=water_dens*g*t_j_depths_unit_iv_bounds(t,2);
        T_t_j_unit_iv_bounds(t,2)=((t_j_depths_unit_iv_bounds(t,2)-intercept)./slope)+273.15;
        D_t_j_unit_iv_bounds(t,2)=4.229e-08.*(exp(-15000./(8.3145.*T_t_j_unit_iv_bounds(t,2))));
        cc_sol_t_j_unit_iv_bounds(t,2)=(sqrt(10.^(-171.9065 - 0.077993.*T_t_j_unit_iv_bounds(t,2) + 2839.319./T_t_j_unit_iv_bounds(t,2) + 71.595.*log10(T_t_j_unit_iv_bounds(t,2)))).*omega_cc).*1000;
        stress_t_j_unit_iv_bounds(t,2)=mean(bulk(depth > unit_boundary_depths(end,1)))*g*t_j_depths_unit_iv_bounds(t,2);
        if use_pff == 1
            sigma_n_t_j_unit_iv_bounds(t,2)=(mean(bulk(depth > unit_boundary_depths(end,1)))-(mean(bulk(depth > unit_boundary_depths(end,1)))*pff))*g*t_j_depths_unit_iv_bounds(t,2);
        else
            sigma_n_t_j_unit_iv_bounds(t,2)=(mean(bulk(depth > unit_boundary_depths(end,1)))-water_dens)*g*t_j_depths_unit_iv_bounds(t,2);
        end
        ps_sr_t_j_unit_iv_bounds(t,2)=(a_d.*((D_t_j_unit_iv_bounds(t,2).*S.*cc_sol_t_j_unit_iv_bounds(t,2))./(d^3)).*((sigma_n_t_j_unit_iv_bounds(t,2).*omega_cc)./(R.*T_t_j_unit_iv_bounds(t,2)))).*(1./((q-(2.*interp1(depth,phi,unit_boundary_depths(end,1)))).^2));
    end
end

if use_rel_vol==2
    for t=1:numel(t_j)
        for z=1:numel(depth)
            ps_sr_t_j(t,z)=ps_sr_t_j(t,z)*cc_rel_vol(z);
        end
        ps_sr_t_j_unit_iv_bounds(t,1)=ps_sr_t_j_unit_iv_bounds(t,1)*interp1(depth,cc_rel_vol,unit_boundary_depths(1,1));
        ps_sr_t_j_unit_iv_bounds(t,2)=ps_sr_t_j_unit_iv_bounds(t,2)*interp1(depth,cc_rel_vol,unit_boundary_depths(end,1));
    end
end

%% integrate with time to get strain

% pre-allocate variables for speed
%first for full depth 
t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
step_t=stepsize*(60*60*24*365.25*1e6);
L_t_j(1:numel(t_j),1:numel(depth))=NaN;
cum_t_j_strain(1:numel(t_j),1:numel(depth))=NaN;
%then for unit IV bounds
t_j_strain_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
L_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
cum_t_j_strain_unit_iv_bounds(1:numel(t_j),1:2)=NaN;
L0=1;

for z=1:numel(depth)
    for t=1:length(t_j)
        t_j_strain(t,z)=-ps_sr_t_j(t,z)*step_t;
        if t==1
            L_t_j(t,z)=L0+(t_j_strain(t,z)*L0);
        elseif isnan(L_t_j(t-1,z))==1
            L_t_j(t,z)=L0+(t_j_strain(t,z)*L0);
        else
            L_t_j(t,z)=L_t_j(t-1,z)+(t_j_strain(t,z)*L_t_j(t-1,z));
        end
        cum_t_j_strain(t,z)=-(L_t_j(t,z)-L0)/L0;
    end
end

final_t_j_strain=cum_t_j_strain(end,:);

for z=1:2
    for t=1:length(t_j)
        t_j_strain_unit_iv_bounds(t,z)=-ps_sr_t_j_unit_iv_bounds(t,z)*step_t;
        if t==1
            L_unit_iv_bounds(t,z)=L0+(t_j_strain_unit_iv_bounds(t,z)*L0);
        elseif isnan(L_unit_iv_bounds(t-1,z))==1
            L_unit_iv_bounds(t,z)=L0+(t_j_strain_unit_iv_bounds(t,z)*L0);
        else
            L_unit_iv_bounds(t,z)=L_unit_iv_bounds(t-1,z)+(t_j_strain_unit_iv_bounds(t,z)*L_unit_iv_bounds(t-1,z));
        end
        cum_t_j_strain_unit_iv_bounds(t,z)=-(L_unit_iv_bounds(t,z)-L0)/L0;
    end
end


for i=1:numel(recovery)
    mean_section_model_strains(i)=mean(final_t_j_strain(depth >= top_cored(i) & depth <= bottom_cored(i)));
end

%% plot conditions and model plots

fd_fontsize = 17;

figure('Units','normalized','OuterPosition',[0.1 0 0.9 1],'Color','w');
subplot(3,2,1)
plot(t_j,t_j_depths_unit_iv_bounds(:,2),'-k',t_j,t_j_depths_unit_iv_bounds(:,1),'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
ylabel('Sediment depth (mbsf)')
legend('Bottom of Unit IV','Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize)
xlim([0 65])
ylim([0 900])
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'FontSize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 900])
set(gca,'box','off')

subplot(3,2,3)
plot(t_j,T_t_j_unit_iv_bounds(:,2)-273.15,'-k',t_j,T_t_j_unit_iv_bounds(:,1)-273.15,'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
ylabel('Temperature (\circC)')
legend('Bottom of Unit IV','Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize)
xlim([0 65])
ylim([0 35])
hold on
yyaxis right
plot(t_j,log10(cc_sol_t_j_unit_iv_bounds(:,2)),'-b','DisplayName','Calcite solubility, bottom of Unit IV','LineWidth',1)
ylabel('Log_{10} solubility (m^3/m^3)')
xlim([0 65])
hold off

subplot(3,2,5)
plot(t_j,stress_t_j_unit_iv_bounds(:,2)./1e6,'-k',...
    t_j,stress_t_j_unit_iv_bounds(:,1)./1e6,'-.k',...
    t_j,sigma_n_t_j_unit_iv_bounds(:,2)./1e6,'-r',...
    t_j,sigma_n_t_j_unit_iv_bounds(:,2)./1e6,'-.r','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
xlabel('Time before present (Myr)')
ylabel('Stress (MPa)')
set(gca,'box','off')
xlim([0 65])
%ylim([0 56])
lgd7=legend('Lithostatic, bottom of Unit IV',...
    'Lithostatic, top of Unit IV','\lambda = 0.95, bottom of Unit IV',...
    '\lambda = 0.95, top of Unit IV',...
    'Location','northwest','Box','off','fontsize',fd_fontsize-1);
%lgd7.FontSize=13;
lgd7.NumColumns = 2;
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 56])
hold off

subplot(3,2,2)
plot(t_j,sigma_n_t_j_unit_iv_bounds(:,2)./1e6,'-r',...
    t_j,sigma_n_t_j_unit_iv_bounds(:,1)./1e6,'-.r','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
%xlabel('Time before present (Myr)')
ylabel({'Effective normal stress','(MPa)'})
xlim([0 65])
leg01=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','west','Box','off','fontsize',fd_fontsize);
%lgd07.FontSize=13;
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([0 9])
set(gca,'box','off')

subplot(3,2,4)
plot(t_j,log10(ps_sr_t_j_unit_iv_bounds(:,2)),'-k',t_j,log10(ps_sr_t_j_unit_iv_bounds(:,1)),'-.k','LineWidth',1)
set(gca,'XDir','reverse','TickDir','out','fontsize',fd_fontsize,'box','off')
%xlabel('Time before present (Myr)')
ylabel('Log_{10} strain rate (s^{-1})')
lgd08=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','northwest','Box','off','fontsize',fd_fontsize);
xlim([0 65])
%lgd08.FontSize=13;
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'fontsize',fd_fontsize,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','YTickLabel',[])
ylim([-20 -16])
set(gca,'box','off')

subplot(3,2,6)
plot(t_j,cum_t_j_strain_unit_iv_bounds(:,2),'-k',t_j,cum_t_j_strain_unit_iv_bounds(:,1),'-.k','LineWidth',1)
set(gca,'XDir','reverse','fontsize',fd_fontsize,'TickDir','out','box','off')
xlabel('Time before present (Myr)')
ylabel('Cumulative strain')
lgd09=legend('Bottom of Unit IV',...
    'Top of Unit IV',...
    'Location','northwest','Box','off',...
    'AutoUpdate','off','fontsize',fd_fontsize);
xlim([0 65])
hold on
%mini plot to add
mini_xmin = 0;
mini_xmax = 1.9;
mini_ymin = 0;
mini_ymax = 9e-4;

plot([mini_xmax mini_xmin mini_xmin mini_xmax mini_xmax],...
    [mini_ymin mini_ymin mini_ymax mini_ymax mini_ymin],'r','LineWidth',1.5)
% plot([mini_xmax 13.55],[0 7e-5],'r','LineWidth',1)
% plot([mini_xmin 1.05],[mini_ymax 5.1e-4],'r','LineWidth',1)
%text(mini_xmax,mini_ymax,'Location of inset','FontSize',15,'HorizontalAlignment','right','VerticalAlignment','bottom')
hold off

%lgd09.FontSize=13;
ax2_pos=get(gca,'Position');
axes('Position',ax2_pos,'YAxisLocation','right','XColor','none','Color','none','TickDir','out','fontsize',fd_fontsize,'YTickLabel',[])
ylim([0 0.015])
set(gca,'box','off')
mini_plot_ax = axes('Position',[0.835 0.125 0.065 0.065],...[ax2_pos(1)*1.035 0.2 0.1 0.12],...
    'Color','none');
plot(t_j,cum_t_j_strain_unit_iv_bounds(:,1),'-.k','LineWidth',1)
xlim([mini_xmin mini_xmax])
ylim([mini_ymin mini_ymax])
mini_plot_ax.YTickLabel{2}=num2str(str2double(mini_plot_ax.YTickLabel{2})*(10^(-4)));
%set(gca,'XDir','reverse','fontsize',15,'XAxisLocation','top',)
mini_plot_ax.XDir = 'reverse';
mini_plot_ax.FontSize =fd_fontsize;
mini_plot_ax.XAxisLocation = 'top';
mini_plot_ax.TickDir = 'out';
mini_plot_ax.XColor = 'red';
mini_plot_ax.YColor = 'red';
mini_plot_ax.LineWidth = 1.2;

tightfig;

if d==1.5e-4
    final_strain_150=final_t_j_strain;
end

if d == saving_size
    export_fig 'conditions_and_model_strain_high_PF.eps' '-painters' 
end

%% check with strains from stylolite
temp=final_t_j_strain;

%get mean strains for small grain size model
mean_unit_model_strain(1,1)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(2,1)));
mean_unit_model_strain(1,2)=mean(temp(depth >= unit_boundary_depths(2,1) & depth < unit_boundary_depths(3,1)));
mean_unit_model_strain(1,3)=mean(temp(depth >= unit_boundary_depths(3,1) & depth < unit_boundary_depths(4,1)));
mean_unit_model_strain(1,4)=mean(temp(depth >= unit_boundary_depths(4,1) & depth < unit_boundary_depths(5,1)));
mean_unit_model_strain(1,5)=mean(temp(depth >= unit_boundary_depths(5,1) & depth < unit_boundary_depths(6,1)));
mean_unit_model_strain(1,6)=mean(temp(depth >= unit_boundary_depths(6,1) & depth < unit_boundary_depths(7,1)));
mean_unit_model_strain(1,7)=mean(temp(depth >= unit_boundary_depths(7,1) & depth < unit_boundary_depths(8,1)));
mean_unit_model_strain(1,8)=mean(temp(depth >= unit_boundary_depths(1,1) & depth < unit_boundary_depths(8,1)));

mean_unit_model_strain(2,:)=mean_unit_model_strain(1,:);

%for plotting out lithological boundaries
xlimss=[10^(-8) 10^(-3.75)];

%% compare strains to those from stylolites
%limits for axes
plot_points=[-4 -1];

figure('units','normalized','outerposition',[1 0 0.8 0.6]);
subplot(1,2,1)
plot(log10(whole_stack_strains),log10(mean_unit_model_strain(1,:)),'.k','MarkerSize',25)
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','out')
legend(['d= ' num2str(d(1)*1e6) ' \mum'])
legend('Box','off','Location','southeast')
if use_rel_vol==1
    title('Comparing XRD-weighted modelled and stylolite based strains')
else
    title('Comparing modelled and stylolite strains')
end
hold on
label_shift=0.08;
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of calculated and modelled strains')
litho_stack_names=litho_names;
litho_stack_names(8)="Total Unit IV";
text(log10(whole_stack_strains)+label_shift,log10(mean_unit_model_strain(1,:)),litho_stack_names,'FontSize',15)%,'Rotation',45)
box off
hold off

for i=1:numel(recovery)
    mean_section_model_strains(i)=mean(final_t_j_strain(depth >= top_cored(i) & depth <= bottom_cored(i)));
end

plot_points=[-4.5 -1];

subplot(1,2,2)
scatter(log10(section_strains),log10(mean_section_model_strains),50,1:1:numel(recovery),'filled')
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','out')
xlim([-4.5 -1.5])
ylim([-4.5 -1.5])
colormap(flip(parula(45)))
c=colorbar('eastoutside');
col_ticks(1)=min(c.Ticks);
col_ticks(2)=max(c.Ticks);
c.Ticks=[1 numel(recovery)];
c.TickLabels={'Top of core','Bottom of core'};
c.Direction='reverse';
c.Label.FontSize=15;
hold on
plot(plot_points,plot_points,'-k')
hold off

%% trigger for grain sizes and plot
if d==1e-5
    litho_strains_10=mean_unit_model_strain(1,:);
    section_strains_10=mean_section_model_strains;
    close all
elseif d==5e-5
    litho_strains_50=mean_unit_model_strain(1,:);
    section_strains_50=mean_section_model_strains;
elseif d==1e-4
    litho_strains_100=mean_unit_model_strain(1,:);
    section_strains_100=mean_section_model_strains;
    close all
end


if d == 1e-4
%% plot out comparison 

plot_points(1)=-4.5;
plot_points(2)=-0.5;

cmp=parula(10);
figure('units','normalized','outerposition',[1 0 0.8 0.6],'Color','w');
subplot(1,2,1)
plot(log10(whole_stack_strains),log10(litho_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(whole_stack_strains),log10(litho_strains_50),'.','Color',cmp(5,:),'MarkerSize',20)
plot(log10(whole_stack_strains),log10(litho_strains_100),'.','Color',cmp(8,:),'MarkerSize',20)
% Highlight total strain value for Unit IV
scatter(log10(whole_stack_strains(end)),log10(litho_strains_10(end)),20,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_50(end)),20,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_100(end)),20,'k','o')
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','in')
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size','Unit IV total')
legend('Box','off','Location','southeast','AutoUpdate','off')
hold on
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains')

litho_strains_plotting(:,1)=rot90(litho_strains_10,3);
litho_strains_plotting(:,2)=rot90(litho_strains_10,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(1,:))
end

litho_strains_plotting(:,1)=rot90(litho_strains_50,3);
litho_strains_plotting(:,2)=rot90(litho_strains_50,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(5,:))
end

litho_strains_plotting(:,1)=rot90(litho_strains_100,3);
litho_strains_plotting(:,2)=rot90(litho_strains_100,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(8,:))
end
hold off

subplot(1,2,2)
plot(log10(section_strains),log10(section_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(section_strains),log10(section_strains_50),'.','Color',cmp(5,:),'MarkerSize',20)
plot(log10(section_strains),log10(section_strains_100),'.','Color',cmp(8,:),'MarkerSize',20)
xlabel('Log_{10} strain from stylolite mass loss and frequency')
ylabel('Log_{10} strain from ps-dissolution model')
set(gca,'fontsize', 15,...
    'TickDir','in')
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size')
legend('Box','off','Location','southeast','AutoUpdate','off')
hold on
label_shift=0.08;
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains')

section_strains_plotting(:,1)=rot90(section_strains_10,3);
section_strains_plotting(:,2)=rot90(section_strains_10,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(1,:))
end

section_strains_plotting(:,1)=rot90(section_strains_50,3);
section_strains_plotting(:,2)=rot90(section_strains_50,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(5,:))
end

section_strains_plotting(:,1)=rot90(section_strains_100,3);
section_strains_plotting(:,2)=rot90(section_strains_100,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(8,:))
end

hold off
tightfig;

%% plot images and final comparison
%create time-depth colour figures
cmp1=parula(20);
fd_plot=figure('Units','normalized','OuterPosition',[0.1 0 0.8 1],'Color','w');
axes('Position',[0.05 0.7 0.35 0.25])
p=pcolor(t_j,depth,transpose(sigma_n_t_j./1e6));
shading interp
colormap(gca,cmp1)
set(p,'AlphaData',transpose(~isnan(sigma_n_t_j)))
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
box off
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
c=colorbar('south');
c.Label.String = 'Effective Normal stress (MPa)';
c.Label.FontSize = fd_fontsize;
c.FontSize=fd_fontsize;
c.TickDirection='out';
c.Position=[0.075 0.92 0.25 0.025];
caxis([0 8])
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

cmp2=parula(20);

axes('Position',[0.455 0.7 0.35 0.25])
p=pcolor(t_j,depth,transpose(log10(cc_sol_t_j)));
shading interp
colormap(gca,cmp2)
set(p,'AlphaData',transpose(~isnan(sigma_n_t_j)))
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
box off
c=colorbar('south');
c.Label.String = 'Log_{10} Calcite solubility (m^3/m^3)';
c.Label.FontSize = fd_fontsize;
c.TickDirection='out';
c.Position=[0.47 0.92 0.25 0.025];
c.FontSize=fd_fontsize;
%c.Direction = 'reverse';
caxis([-5.715 -5.615])
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

cmp3=parula(14);

axes('Position',[0.05 0.06 0.35 0.575])
b=pcolor(t_j,depth,transpose(log10(ps_sr_t_j)));
shading interp
colormap(gca,cmp3)
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
box off
c=colorbar('south');
c.Label.String = 'Log_{10} modelled strain rate (s^{-1})';
c.Label.FontSize = fd_fontsize;
c.Position=[0.075 0.59 0.25 0.025];
c.TickDirection='out';
c.FontSize=fd_fontsize;
if d==5e-5
    caxis([-21 -14.5])
elseif d==1e-4
    caxis([-19 -15.5])
end
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

axes('Position',[0.455 0.06 0.35 0.575])
pre=pcolor(t_j,depth,transpose(cum_t_j_strain));
shading interp
xlabel('Time before present (Myr)')
ylabel('Depth (mbsf)')
set(gca,'YDir','reverse','XDir','reverse','fontsize',fd_fontsize)
ylim([min(unit_boundary_depths(:,1)) max(unit_boundary_depths(:,1))])
xlim([0 65])
if d==5e-5
    caxis([0 0.35])
colormap(gca,parula(14))
elseif d==1e-4
    caxis([0 0.02])
colormap(gca,parula(16))
end
box off
c=colorbar('south');
c.Label.String = 'Modelled cumulative strain';
c.Label.FontSize = fd_fontsize;
c.Position=[0.47 0.59 0.25 0.025];
c.TickDirection='out';
c.FontSize=fd_fontsize;
hold on
plot([max(t_j) interp1(y_short,age_short,unit_boundary_depths(2,1))+3],unit_boundary_depths(2,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(3,1))+8],unit_boundary_depths(3,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(4,1))+7],unit_boundary_depths(4,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(5,1))+2],unit_boundary_depths(5,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(6,1))+2],unit_boundary_depths(6,:),'-k',...
    [max(t_j) interp1(y_short,age_short,unit_boundary_depths(7,1))+2],unit_boundary_depths(7,:),'-k')
hold off

export_fig 'cum_final_condition_images_high_PF.png' -r300
% close

%% calculated concentrated clay per metre at given point from model
clear conc_dist conc_dist_removed clay_conc
conc_dist=1; %distance over which to calculate material removal in m
conc_dist_removed=(rot90(final_t_j_strain,3).*conc_dist);% distance of material removed in m (conc_dist./rot90((1-final_time_fut_strain),3))-conc_dist;
clay_conc=(clay_rel_vol.*conc_dist.*5.5.*5.5)...volume of clay in one metre of core (m3)
    ./((conc_dist-conc_dist_removed).*5.5.*5.5);%volume of remaining material per metre core (m3)
% diff_conc_orig=(clay_conc-clay_rel_vol);
% diff_conc_orig(diff_conc_orig<0)=0;
%get cumulative styl numbers per metre
max_depth=850;
min_depth=501;
cum_dist=1;

cum_styl(:,1:ceil(((max_depth-min_depth)/cum_dist)))=nan;

for i=1:ceil(((max_depth-min_depth)/cum_dist)+1)
    
    cum_styl(i)=sum(styl_depths<=(min_depth+(i*cum_dist))); 
    
end

cum_styl_grad(:,1:numel(cum_styl)-1)=nan;
cum_grad_depth(:,1:numel(cum_styl)-1)=nan;
depths=min_depth:1:max_depth;

for i=1:numel(cum_styl)
    if i<numel(cum_styl)
        cum_styl_grad(i)=cum_styl(i+1)-cum_styl(i);
    else
        cum_styl_grad(i)=0;
    end
    cum_grad_depth(i)=depths(i)+(cum_dist/2)+cum_dist;
    
    
end
cac03_p_m=interp1(depth_carb_unique,carb_unique,(depths-0.5));
strain=0.5;
styl_thickness_p_m=cum_styl_grad.*(styl_width*seams_per_stylolite);
stretch=1-strain;
lo_styl_p_m=styl_thickness_p_m./stretch;
undeformed_thicknesses_p_m=(conc_dist-styl_thickness_p_m)+lo_styl_p_m;
strain_p_m=1-(conc_dist./undeformed_thicknesses_p_m);
    
styl_dist_removed=((strain_p_m)*conc_dist);% distance of material removed in m (conc_dist./rot90((1-final_time_fut_strain),3))-conc_dist;
clay_styl=(interp1(depth,clay_rel_vol,(depths-0.5)).*conc_dist.*5.5.*5.5)...volume of clay in one metre of core (m3)
    ./((conc_dist-styl_dist_removed).*5.5.*5.5);

%% plot out comparison with final strains with depth

fd_fontsize = 17;

plot_points(1)=-5.5;
plot_points(2)=0;
cmp=parula(10);

figure('units','normalized','outerposition',[0 1 0.8 0.6],'Color','w');
subplot(2,5,[1 6])
plot(800,0.06,'Color','none')
hold on
plot(final_strain_150,depth,'-','Color',cmp(8,:),'LineWidth',1.5)
plot(final_t_j_strain,depth,'-','Color',cmp(6,:),'LineWidth',1.5)
plot(strain_p_m,(depths-0.5),'-','Color','k','LineWidth',1)
xlabel('Final strain per metre');

ylabel('Depth (mbsf)');

xticks(0:0.025:0.4)
set(gca,'YDir','reverse',...
    'XAxisLocation','bottom',...
    'fontsize', fd_fontsize,...
    'TickDir','in',...
    'YAxisLocation','left')
ylim([unit_boundary_depths(1,1) unit_boundary_depths(end,1)])
legend('Modelled',[num2str(150) ' \mum grain size'],[num2str(d*1e6) ' \mum grain size'],...
    'Observed','Box','Off','Location','northeast','FontSize',fd_fontsize)


subplot(2,5,[2 8])
plot(log10(whole_stack_strains),log10(litho_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(whole_stack_strains),log10(litho_strains_50),'.','Color',cmp(4,:),'MarkerSize',20)
plot(log10(whole_stack_strains),log10(litho_strains_100),'.','Color',cmp(6,:),'MarkerSize',20)
% Highlight total strain value for Unit IV
scatter(log10(whole_stack_strains(end)),log10(litho_strains_10(end)),40,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_50(end)),40,'k','o')
scatter(log10(whole_stack_strains(end)),log10(litho_strains_100(end)),40,'k','o')
xlabel('Log_{10} strain from stylolite mass loss and frequency');
T = ylabel('Log_{10} strain from pressure solution model');
T.VerticalAlignment = 'baseline';
T.Position(1) = -6.035;
T.Position(2) = -2.25;
%T.Position = ylab_pos;
set(gca,'fontsize', fd_fontsize)
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size','Unit IV total')
legend('Box','off','Location','southeast','AutoUpdate','off','FontSize',fd_fontsize)
hold on
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains','LineWidth',1.5)

xlim(plot_points)
ylim(plot_points)
yticks(-5:1:0)

litho_strains_plotting(:,1)=rot90(litho_strains_10,3);
litho_strains_plotting(:,2)=rot90(litho_strains_10,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(1,:),'LineWidth',1)
end

litho_strains_plotting(:,1)=rot90(litho_strains_50,3);
litho_strains_plotting(:,2)=rot90(litho_strains_50,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(4,:),'LineWidth',1)
end

litho_strains_plotting(:,1)=rot90(litho_strains_100,3);
litho_strains_plotting(:,2)=rot90(litho_strains_100,3);

for i=1:length(plotting_strain_ranges_litho)
plot(log10(plotting_strain_ranges_litho(i,:)),log10(litho_strains_plotting(i,:)),'-','Color',cmp(6,:),'LineWidth',1)
end

line_lab = text(-4.5,-4.5,'Agreement of strains',...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',fd_fontsize);
line_lab.Rotation = 49;
hold off


subplot(2,5,[4 10])
plot(log10(section_strains),log10(section_strains_10),'.','Color',cmp(1,:),'MarkerSize',20)
hold on
plot(log10(section_strains),log10(section_strains_50),'.','Color',cmp(4,:),'MarkerSize',20)
plot(log10(section_strains),log10(section_strains_100),'.','Color',cmp(6,:),'MarkerSize',20)
xlabel('Log_{10} strain from stylolite mass loss and frequency');
set(gca,'fontsize', fd_fontsize)
legend('10 \mum grain size','50 \mum grain size','100 \mum grain size')
legend('Box','off','Location','southeast','AutoUpdate','off','FontSize',fd_fontsize)
hold on
plot(plot_points,plot_points,'-k',...
    'DisplayName','Agreement of strains','LineWidth',1.5)
label_shift=0.08;
line_lab = text(-4.5,-4.5,'Agreement of strains',...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',fd_fontsize);
line_lab.Rotation = 49;
xlim(plot_points)
ylim(plot_points)
yticks(-5:1:0)
section_strains_plotting(:,1)=rot90(section_strains_10,3);
section_strains_plotting(:,2)=rot90(section_strains_10,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(1,:),'LineWidth',1)
end

section_strains_plotting(:,1)=rot90(section_strains_50,3);
section_strains_plotting(:,2)=rot90(section_strains_50,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(4,:),'LineWidth',1)
end

section_strains_plotting(:,1)=rot90(section_strains_100,3);
section_strains_plotting(:,2)=rot90(section_strains_100,3);

for i=1:length(plotting_strain_ranges_section)
plot(log10(plotting_strain_ranges_section(i,:)),log10(section_strains_plotting(i,:)),'-','Color',cmp(6,:),'LineWidth',1)
end


hold off
tightfig;
export_fig 'final_strain_comparison_var_grain_size_high_PF.eps' '-painters' 
end
end
% close
