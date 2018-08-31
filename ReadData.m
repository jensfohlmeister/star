
function[depth, fm, fm_err, Anchorpoint_age, Anchorpoint_depth] = ReadData(key)

% reads data or produces a random dataset


f=figure(100);
set(f, 'Units','centimeters', 'Position', [3 12 20 13])
hold on

%%% reading data of your data set
if key == 5 
    
    [value] = textread('parameters.txt',...
        '%*s %*s %s %*s', 'commentstyle','matlab');%

    Anchorpoint_age     = str2num(value{2});
    Anchorpoint_depth   = str2num(value{3}); 

    stalfile    = value{6}; 
    proxyfile_d13C   = value{7};
    proxyfile_MgCa   = value{8};
    
    [depth, fm, fm_err] = textread(stalfile,...
    '%f %f %f','commentstyle','matlab'); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% check for trends in ancillary proxy data (d13C, Mg/Ca, etc.) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Will only be run if the proxy data is available, otherwise the program
    %continues without this information
    if strcmpi(proxyfile_d13C,'none') == 1,
        disp('No ancillary d13C proxy data available')
    else
        Trends_d13C(proxyfile_d13C);
    end

    if strcmpi(proxyfile_MgCa,'none') == 1,
        disp('No ancillary MgCa proxy data available')
    else
        Trends_MgCa(proxyfile_MgCa);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

%%% our artificial data set (see paper);
elseif key == 4
    [depth, fm, fm_err] = textread('dataArtificial.txt',...
    '%f %f %f','commentstyle','matlab');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting the age depth model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_depth = [0 0.1 3 7.5 12.5 17.5 22.5 27.5 32.5 37.5 42.5 47.5 57.5 67.5 77.5 87.5 97.5 107.5 117.5 127.5 137.5 147.5 157.5 167.5 177.5 192.5]; 
    orig_age = [98 100 147 221 302 384 465 547 629 710 792 873 1036 1200 1363 2115 2866 3618 4370 5122 5874 6626 7378 8130 8881 10009];

    figure(2)
    hold on
    box on
    plot(orig_depth, orig_age, 'o-k','LineWidth',2)
    plot([orig_depth(15) orig_depth(15)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    ylabel('age [a]','Fontsize',14,'FontWeight','bold')
    xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
    set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define anchor point (and uncertainties?) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Anchorpoint_age = [98 105 90]; % [Age UpperError LowerError]
    Anchorpoint_depth = 0; %[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% randomly constructed, artificial dataset, 1 breakbpoint in growth rate
elseif key == 1


    stal_length = round(400 + (600-400)*rand); % in [mm]
    number_14C_Datapoints = round( (5+(10-5)*rand) *stal_length/100); % per 100mm of stalagmite 5-10 14C measurements

    msd = stal_length/number_14C_Datapoints; %mean sampling distance
    depth(:,1) = [0.5:msd:stal_length stal_length]; % depths of 14C sampling
    depth(2:end-1,1) = round( (depth(2:end-1,1) + (-msd/2 + msd*rand(length(depth)-2,1) ))*10 )/10; %make depth it a bit more random

    %%% defining two mean growth rates and breakpoint (make it all to 'rand'later on)
    GR1 = 0.0923; %[mm/a]
    GR2 = 0.0621; %[mm/a]
    bp = stal_length/3+ (2/3*stal_length-stal_length/3)*rand; % remove the boundaries after testing period!!!
    index_bp = find( depth<bp,1,'last' );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define anchor point (and uncertainties?) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Anchorpoint_age = [100 100 100]; % [Age UpperError LowerError]
    Anchorpoint_depth = depth(1); %[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% find age for all depths

    age(1,1)=Anchorpoint_age(1); %[years]
    for i = 2:index_bp
        age(i,1) = age(1,1) + (depth(i,1)-depth(1,1))/GR1;
    end
    age_bp = age(1,1) + (bp-depth(1,1))/GR1;
    for i = index_bp+1:length(depth)
        age(i,1) = age_bp + (depth(i,1)-bp)/GR2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting the age depth model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_depth = [depth(1:index_bp)' bp depth(index_bp+1:length(depth))'];
    orig_age = [age(1:index_bp)' age_bp age(index_bp+1:length(depth))'];

    figure(2)
    hold on
    box on
    plot(orig_depth, orig_age, 'k','LineWidth',2)
    plot(depth, age, 'ok')
    plot([bp bp], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    ylabel('age [a]','Fontsize',14,'FontWeight','bold')
    xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
    set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% read IntCal13, for determination of atm. a14C values at defined
    %%% ages and calculation of stalagmite 14C value
    [Intcal_age, c14age, errC14age, Intcal_D14C, ~] = textread('intcal13.txt',...
    '%f %f %f %f %f','commentstyle','matlab');
    for i = 1:length(age)
        index_fm = find(abs(Intcal_age-age(i)) == min(abs( Intcal_age-age(i))) );
        fm_atm(i,1) = exp(-c14age(index_fm)/8033);
        DCF(i,1) = 0.12 + (0.17-0.12)*rand; %define random DCF between 0.12 and 0.18
        fm(i,1) = (1-DCF(i))*fm_atm(i); % artificial stalagmite 14C
        fm_err(i,1) = fm(i)*0.003; %three permil error for all samples
    end

%%% randomly constructed, artificial dataset, 2 breakbpoint in growth rate   
elseif key == 2 

    stal_length = round(100 + (800-100)*rand); % in [mm]
    number_14C_Datapoints = round( (5+(10-5)*rand) *stal_length/100); % per 100mm of stalagmite 5-10 14C measurements

    msd = stal_length/number_14C_Datapoints; %mean sampling distance
    depth(:,1) = [0.5:msd:stal_length stal_length]; % depths of 14C sampling
    depth(2:end-1,1) = round( (depth(2:end-1,1) + (-msd/2 + msd*rand(length(depth)-2,1) ))*10 )/10; %make depth it a bit more random

    %%% defining two mean growth rates and breakpoint (make it all to 'rand'later on)
    GR1 = 0.1221; %[mm/a]
    GR2 = 0.0554; %[mm/a]
    GR3 = 0.221; %[mm/a]
    bp(1) = stal_length/5+ (2/5*stal_length-stal_length/5)*rand; % remove the boundaries after testing period!!!
    bp(2) = stal_length*3/5+ (4/5*stal_length-stal_length*3/5)*rand; % remove the boundaries after testing period!!!
    
    index_bp(1) = find( depth<bp(1),1,'last' );
    index_bp(2) = find( depth<bp(2),1,'last' );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define anchor point (and uncertainties?) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Anchorpoint_age = [100 100 100]; % [Age UpperError LowerError]
    Anchorpoint_depth = depth(1); %[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% find age for all depths

    age(1,1)=Anchorpoint_age(1); %[years]
    for i = 2:index_bp(1)
        age(i,1) = age(1,1) + (depth(i,1)-depth(1,1))/GR1;
    end
    age_bp(1) = age(1,1) + (bp(1)-depth(1,1))/GR1;
    for i = index_bp(1)+1:index_bp(2)
        age(i,1) = age_bp(1) + (depth(i,1)-bp(1))/GR2;
    end
    age_bp(2) = age_bp(1) + (bp(2)-bp(1))/GR2;
    for i = index_bp(2)+1:length(depth)
        age(i,1) = age_bp(2) + (depth(i,1)-bp(2))/GR3;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting the age depth model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_depth = [depth(1:index_bp(1))' bp(1) depth(index_bp(1)+1:index_bp(2))' bp(2) depth(index_bp(2)+1:length(depth))'];
    orig_age = [age(1:index_bp(1))' age_bp(1) age(index_bp(1)+1:index_bp(2))' age_bp(2) age(index_bp(2)+1:length(depth))'];

    figure(2)
    hold on
    box on
    plot(orig_depth, orig_age, 'k','LineWidth',2)
    plot(depth, age, 'ok')
    plot([bp(1) bp(1)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    plot([bp(2) bp(2)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    ylabel('age [a]','Fontsize',14,'FontWeight','bold')
    xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
    set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% read IntCal13, for determination of atm. a14C values at defined
    %%% ages and calculation of stalagmite 14C value
    [Intcal_age, c14age, errC14age, Intcal_D14C, ~] = textread('intcal13.txt',...
    '%f %f %f %f %f','commentstyle','matlab');
    for i = 1:length(age)
        index_fm = find(abs(Intcal_age-age(i)) == min(abs( Intcal_age-age(i))) );
        fm_atm(i,1) = exp(-c14age(index_fm)/8033);
        DCF(i,1) = 0.12 + (0.18-0.12)*rand; %define random DCF between 0.12 and 0.18
        fm(i,1) = (1-DCF(i))*fm_atm(i); % artificial stalagmite 14C
        fm_err(i,1) = fm(i)*0.003; %three permil error for all samples
    end
    
%%% randomly constructed, artificial dataset, 3 breakbpoint in growth rate
elseif key == 3 

    stal_length = round(100 + (800-100)*rand); % in [mm]
    number_14C_Datapoints = round( (5+(10-5)*rand) *stal_length/100); % per 100mm of stalagmite 5-10 14C measurements

    msd = stal_length/number_14C_Datapoints; %mean sampling distance
    depth(:,1) = [0.5:msd:stal_length stal_length]; % depths of 14C sampling
    depth(2:end-1,1) = round( (depth(2:end-1,1) + (-msd/2 + msd*rand(length(depth)-2,1) ))*10 )/10; %make depth a bit more random

    %%% defining two mean growth rates and breakpoint (make it all to 'rand'later on)
    GR1 = 0.0721; %[mm/a]
    GR2 = 0.254; %[mm/a]
    GR3 = 0.071; %[mm/a]
    GR4 = 0.2; %[mm/a]
    bp(1) = stal_length/7+ (2/7*stal_length-stal_length/7)*rand; % remove the boundaries after testing period!!!
    bp(2) = stal_length*3/7+ (3/7*stal_length-stal_length*2/7)*rand; % remove the boundaries after testing period!!!
    bp(3) = stal_length*5/7+ (6/7*stal_length-stal_length*5/7)*rand; % remove the boundaries after testing period!!!

    index_bp(1) = find( depth<bp(1),1,'last' );
    index_bp(2) = find( depth<bp(2),1,'last' );
    index_bp(3) = find( depth<bp(3),1,'last' );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define anchor point (and uncertainties?) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Anchorpoint_age = [100 100 100]; % [Age UpperError LowerError]
    Anchorpoint_depth = depth(1); %[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% find age for all depths

    age(1,1)=Anchorpoint_age(1); %[years]
    for i = 2:index_bp(1)
        age(i,1) = age(1,1) + (depth(i,1)-depth(1,1))/GR1;
    end
    age_bp(1) = age(1,1) + (bp(1)-depth(1,1))/GR1;
    for i = index_bp(1)+1:index_bp(2)
        age(i,1) = age_bp(1) + (depth(i,1)-bp(1))/GR2;
    end
    age_bp(2) = age_bp(1) + (bp(2)-bp(1))/GR2;
    for i = index_bp(2)+1:index_bp(3)
        age(i,1) = age_bp(2) + (depth(i,1)-bp(2))/GR3;
    end
    age_bp(3) = age_bp(2) + (bp(3)-bp(2))/GR3;
    for i = index_bp(3)+1:length(depth)
        age(i,1) = age_bp(3) + (depth(i,1)-bp(3))/GR4;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting the age depth model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_depth = [depth(1:index_bp(1))' bp(1) depth(index_bp(1)+1:index_bp(2))' bp(2) depth(index_bp(2)+1:index_bp(3))' bp(3) depth(index_bp(3)+1:length(depth))'];
    orig_age = [age(1:index_bp(1))' age_bp(1) age(index_bp(1)+1:index_bp(2))' age_bp(2) age(index_bp(2)+1:index_bp(3))' age_bp(3) age(index_bp(3)+1:length(depth))'];
    %Anchorpoint_age = [orig_age(index_bp(2)+4) orig_age(index_bp(2)+4)+100 orig_age(index_bp(2)+4)-100]; % [Age UpperError LowerError]
    %Anchorpoint_depth = orig_depth(index_bp(2)+4); %[mm]

    figure(2)
    hold on
    box on
    plot(orig_depth, orig_age, 'k','LineWidth',2)
    plot(depth, age, 'ok')
    plot([bp(1) bp(1)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    plot([bp(2) bp(2)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    plot([bp(3) bp(3)], [orig_age(1) orig_age(end)],'k','LineWidth',2)
    ylabel('age [a]','Fontsize',14,'FontWeight','bold')
    xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
    set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% read IntCal13, for determination of atm. a14C values at defined
    %%% ages and calculation of stalagmite 14C value
    [Intcal_age, c14age, errC14age, Intcal_D14C, ~] = textread('intcal13.txt',...
    '%f %f %f %f %f','commentstyle','matlab');
    for i = 1:length(age)
        index_fm = find(abs(Intcal_age-age(i)) == min(abs( Intcal_age-age(i))) );
        fm_atm(i,1) = exp(-c14age(index_fm)/8033);
        DCF(i,1) = 0.12 + (0.18-0.12)*rand; %define random DCF between 0.12 and 0.18
        fm(i,1) = (1-DCF(i))*fm_atm(i); % artificial stalagmite 14C
        fm_err(i,1) = fm(i)*0.003; %three permil error for all samples
    end
    

end %%% end reading stalagmite 14C data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%Check for trends in ancillary proxy data (Mg/Ca)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fit_MgCa, slope_MgCa]= Trends_MgCa(proxyfile_MgCa)

%%% read stalagmite proxy data and check for suitability of dataset %%%%%%%%%
[depth2, MgCa] = textread(proxyfile_MgCa,...
    '%f %f','commentstyle','matlab');
%Define MgCa vs. depth trend
fit_MgCa = polyfit(depth2, MgCa,1);
slope_MgCa = fit_MgCa(:,1);
gerade_MgCa = fit_MgCa(1)*depth2+fit_MgCa(2);

%If the slope is larger than 0.1 or smaller than -0.1, a warning comes up (a significant trend
%might be present)

if abs(slope_MgCa) > 0.1
    disp('Trend detected in ancillary MgCa proxy data: DCF stationarity might be compromised')
    %Plot MgCa
    figure(11)
    plot(depth2, MgCa, '-r')
    hold on
    plot(depth2, gerade_MgCa, '-b')
end

figure(100)
subplot(2,5,[9,10])
yyaxis left
plot(depth2,MgCa,'-k')
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold','YColor', [0 0 0])
xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
ylabel('Mg/Ca','Fontsize',14,'FontWeight','bold', 'color', [0 0 0])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%Check for trends in ancillary proxy data (d13C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fit_d13C, slope_d13C]= Trends_d13C(proxyfile_d13C)

%%% read stalagmite proxy data and check for suitability of dataset %%%%%%%%%
[depth2, d13C] = textread(proxyfile_d13C,...
    '%f %f','commentstyle','matlab');

%Define d13C vs. depth trend
fit_d13C = polyfit(depth2, d13C,1);
slope_d13C = fit_d13C(:,1);
gerade_d13C = fit_d13C(1)*depth2+fit_d13C(2);

%If the slope is larger than 0.1 or smaller than -0.1, a warning comes up (a significant trend
%might be present)

if abs(slope_d13C) > 0.1
    disp('Trend detected in ancillary d13C proxy data: DCF stationarity might be compromised')
    %Plot d13C
    figure(10)
    hold on
    plot(depth2, d13C, '-r')
    plot(depth2, gerade_d13C, '-b')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(100);
subplot(2,5,[9,10]);
yyaxis right
plot(depth2,d13C, 'color', [0.7 0.7 0.7])
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor', [0.7 0.7 0.7])
xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
ylabel(['\delta^{13}C [',char(8240),']'],'Fontsize',14,'FontWeight','bold',...
    'Interpreter','tex', 'color', [0.7 0.7 0.7] )



