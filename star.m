% STAlagmite dating by Radiocarbon ; [star]     version 02.02.2018
%The age model calculates accurate chronologies for stalagmites using 14C.
%Suitable samples are characterized by: 
%   i) the absence of long-term secular variability in DCF; 
%   ii) long enough growth for measurable 14C decay to have taken place 
%   iv) an anchor point of known age is available.
%
%The input file should contain the following columns :
%Depth from top, measured fm, error in fm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STAlagmite dating by Radiocarbon (star): a software tool for reliable and fast age depth modelling
%
%The software comes along with the following publication: Fohlmeister, J., Lechleitner, F.A., 
%STAlagmite dating by Radiocarbon (star): a software tool for reliable and fast age depth 
%modelling, Quaternary Geochronology, submitted
%
%based on previous work, presented in:
%
%Lechleitner, F. A., Fohlmeister, J., McIntyre, C., Baldini, L., Jamieson, R. A., Hercman, H., 
%Gasiorowski, M., Pawlak, J., Stefaniak, K., Socha, P., Eglinton, T. I., Baldini, J. U. L., 2016. 
%A novel approach for reconstruction of accurate radiocarbon-based chronologies for speleothems. 
%Quaternary Geochronology 35, 54-66.
%
%Please cite the paper(s) if you are using this fantastic software. ;)
%
%Please read the paper before considering this software to construct an age-depth model for your climate archive.
%
%A one-page manual is provided for a quick start. Please read this. Afterwards you should be able 
%to handle the most important features of 'star' and application of your dataset shouldn't be a problem anymore.
%
%Star is open-source software; you are free to use, copy, distribute and modify it. However, 
%it does not come with any warranty and the authors do not assume any responsibility for the 
%usefulness of any portion of this program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

[value1, value2, value, value4] = textread('parameters.txt',...
'%s %s %s %s', 'commentstyle','matlab');%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% read stalagmite 14C data (or build an artificial data set) and 
%%% check for suitability of dataset 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
key                     = str2num(value{1});
NRCV                    = str2num(value{4}); %Number of significant radiocarbon values needed for fitting (growth rate changes)
key_breakpoint_age_diff = str2num(value{5}); %Number of years needed for fitting (growth rate changes);
proxyfile_d13C          = value{7};
proxyfile_MgCa          = value{8};

    %key =  1 random, one breakpoint artificial dataset;
    %       2 random, two breakpoint artificial dataset; 
    %       3 random, three breakpoint artificial dataset; 
    %       4 our artificial data set (see Lechleitner et al., 2016);     
    %       5 your own data set
    
[depth, fm, fm_err, Anchorpoint_age, Anchorpoint_depth] = ReadData(key);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check for depth inversions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check1 = find(diff(depth)<=0);

if isempty(check1)==1
    clear check1;
elseif length(check1) == length(index)
    if check1 == index, clear check1;end 
end
A = exist('check1','var');

if A==1
    str = sprintf('The current model version can only deal with sorted data. \nPlease provide the data sorted by depth from top to bottom.');
    error(str)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Identify and specify hiatuses; 
%%% Growth stop positions are marked and age model is split 
%%% check in which part with respect to the growth stop(s) 
%%% the anchor point lies (Anchorpoint_depth_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[index, Anchorpoint_depth_index] = Hiatus_detection1(fm, depth, Anchorpoint_depth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot of log(fm/A0) vs. depth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% plot(depth(1:index(2)),log(fm(1:index(2))),'o-');
% hold on
% for i = 2:length(index)-1
%         plot(depth(index(i)+1:index(i+1)),log(fm(index(i)+1:index(i+1))),'o-');
% end
% set(gca,'Ydir','reverse');
% title('Fm 14C vs. depth')
% xlabel('Depth from top (mm)')
% ylabel('log(fm/A0)') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fit the growth rate function 
%%% of the stalagmite part including the anchor point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for  i = 1: length(index)
    if Anchorpoint_depth_index == i, break; end
end

if i==1  %if Anchor point in the youngest growth section

elseif i == 2 % if anchor point in the second growth section --> reorder the measurements!
    depth_dummy = depth(index(2)+1:index(3));
    depth_dummy(end+1:end+length(depth(index(1):index(2))))=depth(index(1):index(2));
    fm_dummy = fm(index(2)+1:index(3));
    fm_dummy(end+1:end+length(depth(index(1):index(2))))=fm(index(1):index(2));
    index(2) = index(3)-index(2);
    
    depth(1:index(3)) = depth_dummy;
    fm(1:index(3)) = fm_dummy;
    
    clear fm_dummy depth_dummy

elseif i == 3 % if anchor point in the third growth section --> reorder the measurements!
    depth_dummy = depth(index(3)+1:index(4));
    depth_dummy(end+1:end+length(depth(index(2)+1:index(3))))=depth(index(2)+1:index(3));
    depth_dummy(end+1:end+length(depth(index(1):index(2))))=depth(index(1):index(2));
    fm_dummy = fm(index(3)+1:index(4));
    fm_dummy(end+1:end+length(depth(index(2)+1:index(3))))=fm(index(2)+1:index(3));
    fm_dummy(end+1:end+length(depth(index(1):index(2))))=fm(index(1):index(2));
    index_dummy = index(3) - index(2);
    index(2) = index(4) - index(3);
    index(3) = index(2) + index_dummy;
    
    clear depth fm
    
    depth = depth_dummy;
    fm = fm_dummy;
    
    clear fm_dummy depth_dummy index_dummy

end

i=1;
[log_fm1, gerade1, Intcal_A0, Age_const, Age_lower_ci, Age_upper_ci, lambda]  = ...
    GR_fit(depth, fm, [index(i) index(i+1)+1], Anchorpoint_age, Anchorpoint_depth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking for longer-term DCF trends to find obvious growth rate changes
%%% if enough data points avalaibly and growth section long enough (>500a)
%%% otherwise DCF variability might hide any decay and the fitting result 
%%% is too strongly affected 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(depth(index(1):index(2)))>=2*NRCV & Age_const(end)-Age_const(1)>key_breakpoint_age_diff
        breakpoint = SlopeBreakpoint(log_fm1,depth(index(1):index(2)),gerade1,[index(i) index(i+1)+1],NRCV);
        number_bp = length(breakpoint);
    else
        breakpoint = [ ];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if no longer-term residual 14C trends could be observed (or too few 14C
%%% measurements available), do nothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(breakpoint)==1   
    
    else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% otherwise split the record into 'number_bp' growth sections and repeat 
%%% the fitting procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        hold on
        for ii = 1 : number_bp
            plot([breakpoint(ii) breakpoint(ii)], [Age_const(1) Age_const(index(2)-1)],'g')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            index_bp1(ii) = find(abs(depth-breakpoint(ii))==min(abs(depth-breakpoint(ii))) );
        end

        if Anchorpoint_depth <= depth(index_bp1(1))
        %%% fitting the first part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [1 index_bp1(1)+1], Anchorpoint_age, Anchorpoint_depth);

            Intcal_A0 = Intcal_A0_dummy;
            Age_const = Age_const_dummy;
            Age_lower_ci = Age_lower_ci_dummy;
            Age_upper_ci = Age_upper_ci_dummy;
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy
        
        %%% fitting the middle part
            for ii=2:number_bp
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(end) Age_upper_ci(end) Age_lower_ci(end)], depth(index_bp1(ii-1)));

                Intcal_A0(end+1:end+length(Age_const_dummy)-1) = Intcal_A0_dummy(2:end);
                Age_const(end+1:end+length(Age_const_dummy)-1) = Age_const_dummy(2:end);
                Age_lower_ci(end+1:end+length(Age_const_dummy)-1) = Age_lower_ci_dummy(2:end);
                Age_upper_ci(end+1:end+length(Age_const_dummy)-1) = Age_upper_ci_dummy(2:end);
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            
        %%% fitting the oldest part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index_bp1(end) index(i+1)+1], [Age_const(end) Age_upper_ci(end) Age_lower_ci(end)], depth(index_bp1(end)));
    
            Intcal_A0(end+1:end+length(Age_const_dummy)-1) = Intcal_A0_dummy(2:end);
            Age_const(end+1:end+length(Age_const_dummy)-1) = Age_const_dummy(2:end);
            Age_lower_ci(end+1:end+length(Age_const_dummy)-1) = Age_lower_ci_dummy(2:end);
            Age_upper_ci(end+1:end+length(Age_const_dummy)-1) = Age_upper_ci_dummy(2:end);



        elseif Anchorpoint_depth >= depth(index_bp1(end))
        %%% fitting the oldest part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index_bp1(end) index(i+1)+1], Anchorpoint_age, Anchorpoint_depth);
    
            Intcal_A0 = Intcal_A0_dummy;
            Age_const = Age_const_dummy;
            Age_lower_ci = Age_lower_ci_dummy;
            Age_upper_ci = Age_upper_ci_dummy;
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

        %%% fitting the middle part
            for ii=number_bp:-1:2
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(1) Age_upper_ci(1) Age_lower_ci(1)], depth(index_bp1(ii)));

                Intcal_A0_dummy(end+1:end+length(Age_const)-1) = Intcal_A0(2:end);
                Age_const_dummy(end+1:end+length(Age_const)-1) = Age_const(2:end);
                Age_lower_ci_dummy(end+1:end+length(Age_const)-1) = Age_lower_ci(2:end);
                Age_upper_ci_dummy(end+1:end+length(Age_const)-1) = Age_upper_ci(2:end);
                clear Intcal_A0 Age_const Age_lower_ci Age_upper_ci
                Intcal_A0 = Intcal_A0_dummy;
                Age_const = Age_const_dummy;
                Age_lower_ci = Age_lower_ci_dummy;
                Age_upper_ci = Age_upper_ci_dummy;
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            
        %%% fitting the first part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [1 index_bp1(1)+1], [Age_const(1) Age_upper_ci(1) Age_lower_ci(1)], depth(index_bp1(1)));

            Intcal_A0_dummy(end+1:end+length(Age_const)-1) = Intcal_A0(2:end);
            Age_const_dummy(end+1:end+length(Age_const)-1) = Age_const(2:end);
            Age_lower_ci_dummy(end+1:end+length(Age_const)-1) = Age_lower_ci(2:end);
            Age_upper_ci_dummy(end+1:end+length(Age_const)-1) = Age_upper_ci(2:end);
            clear Intcal_A0 Age_const Age_lower_ci Age_upper_ci
            Intcal_A0 = Intcal_A0_dummy;
            Age_const = Age_const_dummy;
            Age_lower_ci = Age_lower_ci_dummy;
            Age_upper_ci = Age_upper_ci_dummy;
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy
         
        elseif Anchorpoint_depth >= depth(index_bp1(1)) & Anchorpoint_depth <= depth(index_bp1(2))
            
        %%% fitting the middle part
            for ii=2:2
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], Anchorpoint_age, Anchorpoint_depth );

                Intcal_A0 = Intcal_A0_dummy;
                Age_const = Age_const_dummy;
                Age_lower_ci = Age_lower_ci_dummy;
                Age_upper_ci = Age_upper_ci_dummy;
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            for ii=number_bp:-1:3
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(end) Age_upper_ci(end) Age_lower_ci(end)], depth(index_bp1(ii-1)) );

                Intcal_A0(end+1:end+length(Age_const_dummy)-1) = Intcal_A0_dummy(2:end);
                Age_const(end+1:end+length(Age_const_dummy)-1) = Age_const_dummy(2:end);
                Age_lower_ci(end+1:end+length(Age_const_dummy)-1) = Age_lower_ci_dummy(2:end);
                Age_upper_ci(end+1:end+length(Age_const_dummy)-1) = Age_upper_ci_dummy(2:end);
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
           
        %%% fitting the oldest part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index_bp1(end) index(i+1)+1], [Age_const(end) Age_upper_ci(end) Age_lower_ci(end)], depth(index_bp1(end)));
    
            Intcal_A0(end+1:end+length(Age_const_dummy)-1) = Intcal_A0_dummy(2:end);
            Age_const(end+1:end+length(Age_const_dummy)-1) = Age_const_dummy(2:end);
            Age_lower_ci(end+1:end+length(Age_const_dummy)-1) = Age_lower_ci_dummy(2:end);
            Age_upper_ci(end+1:end+length(Age_const_dummy)-1) = Age_upper_ci_dummy(2:end);
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy


         %%% fitting the first part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [1 index_bp1(1)+1], [Age_const(1) Age_upper_ci(1) Age_lower_ci(1)], depth(index_bp1(1)));

            Intcal_A0_dummy(end+1:end+length(Age_const)-1) = Intcal_A0(2:end);
            Age_const_dummy(end+1:end+length(Age_const)-1) = Age_const(2:end);
            Age_lower_ci_dummy(end+1:end+length(Age_const)-1) = Age_lower_ci(2:end);
            Age_upper_ci_dummy(end+1:end+length(Age_const)-1) = Age_upper_ci(2:end);
            clear Intcal_A0 Age_const Age_lower_ci Age_upper_ci
            Intcal_A0 = Intcal_A0_dummy;
            Age_const = Age_const_dummy;
            Age_lower_ci = Age_lower_ci_dummy;
            Age_upper_ci = Age_upper_ci_dummy;
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy
           
        else 
        %%% fitting the middle part
            for ii=number_bp:-1:3
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], Anchorpoint_age, Anchorpoint_depth );
              
                Intcal_A0 = Intcal_A0_dummy;
                Age_const = Age_const_dummy;
                Age_lower_ci = Age_lower_ci_dummy;
                Age_upper_ci = Age_upper_ci_dummy;
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end

        %%% fitting the oldest part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index_bp1(end) index(i+1)+1], [Age_const(end) Age_upper_ci(end) Age_lower_ci(end)], depth(index_bp1(end)));
    
            Intcal_A0(end+1:end+length(Age_const_dummy)-1) = Intcal_A0_dummy(2:end);
            Age_const(end+1:end+length(Age_const_dummy)-1) = Age_const_dummy(2:end);
            Age_lower_ci(end+1:end+length(Age_const_dummy)-1) = Age_lower_ci_dummy(2:end);
            Age_upper_ci(end+1:end+length(Age_const_dummy)-1) = Age_upper_ci_dummy(2:end);
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy


        %%% fitting the middle part
            for ii=2:2
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(1) Age_upper_ci(1) Age_lower_ci(1)], depth(index_bp1(ii)));
 
                Intcal_A0_dummy(end+1:end+length(Age_const)-1) = Intcal_A0(2:end);
                Age_const_dummy(end+1:end+length(Age_const)-1) = Age_const(2:end);
                Age_lower_ci_dummy(end+1:end+length(Age_const)-1) = Age_lower_ci(2:end);
                Age_upper_ci_dummy(end+1:end+length(Age_const)-1) = Age_upper_ci(2:end);
                clear Intcal_A0 Age_const Age_lower_ci Age_upper_ci
                Intcal_A0 = Intcal_A0_dummy;
                Age_const = Age_const_dummy;
                Age_lower_ci = Age_lower_ci_dummy;
                Age_upper_ci = Age_upper_ci_dummy;
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            
         %%% fitting the first part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [1 index_bp1(1)+1], [Age_const(1) Age_upper_ci(1) Age_lower_ci(1)], depth(index_bp1(1)));

            Intcal_A0_dummy(end+1:end+length(Age_const)-1) = Intcal_A0(2:end);
            Age_const_dummy(end+1:end+length(Age_const)-1) = Age_const(2:end);
            Age_lower_ci_dummy(end+1:end+length(Age_const)-1) = Age_lower_ci(2:end);
            Age_upper_ci_dummy(end+1:end+length(Age_const)-1) = Age_upper_ci(2:end);
            clear Intcal_A0 Age_const Age_lower_ci Age_upper_ci
            Intcal_A0 = Intcal_A0_dummy;
            Age_const = Age_const_dummy;
            Age_lower_ci = Age_lower_ci_dummy;
            Age_upper_ci = Age_upper_ci_dummy;
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy
  
            
        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

clear index_bp1    
    


%%% Age determination in sections without anchor point %%%%%%%%%%%%%%%%%%%

for i = 2:length(index)-1
     
     
    %Calculate DCF for the top part % with my procedure and include errors!
    stal_ainit = fm(index(i-1)+1:index(i)-1) ./ exp(-lambda*Age_const(index(i-1)+1:index(i)-1));
    DCF = 1-(stal_ainit./Intcal_A0(index(i-1)+1:index(i)-1));% DCF of the top samples
    avg_DCF_top = mean(DCF); %use this DCF to estimate the "starting point" of the second part (below the hiatus)
    DCF_STD = std(DCF); %Assume constant DCF over the hiatus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% calculate atm. 14C signal at age, determined in previous stalagmite section %%%%%%%%%%%%
    C14_age(index(i-1)+1:index(i)-1) = -8033 * log(fm(index(i-1)+1:index(i)-1));

    [atm14C,atm14C_1sig] = montecarlo(Age_const((index(i-1)+1:index(i)-1))/1000,...
        (Age_upper_ci(index(i-1)+1:index(i)-1)-Age_lower_ci(index(i-1)+1:index(i)-1))/2/2/1000); %1 sigma errors
    
    res_age(index(i-1)+1:index(i)-1) = C14_age(index(i-1)+1:index(i)-1)- atm14C(1:end); %reservoir age
    avg_res_age = mean(res_age(index(i-1)+1:index(i)-1));
    STD_res_age = std(res_age(index(i-1)+1:index(i)-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Apply this Reservoir age in this stalagmite growth stection %%%%%%%
    C14_age(index(i)+1:index(i+1)) = -8033 * log(fm(index(i)+1:index(i+1)));
    Age_corr = C14_age(index(i)+1:index(i+1)) - avg_res_age;
    for k = 1:length(Age_corr)
        C14_age_err(index(i)+k) = abs(-8033 * fm_err(index(i)+k)/fm(index(i)+k));
        Age_corr_err(k) = sqrt(STD_res_age^2 + C14_age_err(index(i)+k)^2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Calibrate DCF corrected 14C values for fitting purposes of the 
    %%% growth rate for the next (older) stalagmite section after the
    %%% hiatus
    fprintf('Average DCF from growth section %g: %g +/- %G %% \n\n', i-1,avg_DCF_top*100,DCF_STD*100)
    fprintf('Average Res. Age from growth section %g: %g +/- %g years \n\n', i-1,avg_res_age,STD_res_age)
    disp('Calibrating the DCF corrected 14C values of the older stalagmite section ...')
    [x2, y2] = calib2(Age_corr, Age_corr_err, depth(index(i)+1:index(i+1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    disp(' ')
    fprintf('Calculating the age-depth model for growth section %g of %g sections \n\n', i,(length(index)-1))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fitting GR in sections without anchor point
    %%% function slightly different compared to GR_fit.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [log_fm1, gerade1, Intcal_A01, Age_const1, Age_lower_ci1, Age_upper_ci1, lambda]  = ...
        GR2_fit(depth, fm, [index(i) index(i+1)+1], Age_const(index(i)), depth(index(i)), x2, y2);

    Intcal_A0(index(i)+1:index(i+1)) = Intcal_A01(index(i)+1:index(i+1));
    Age_const(index(i)+1:index(i+1)) = Age_const1(index(i)+1:index(i+1));
    Age_lower_ci(index(i)+1:index(i+1)) = Age_lower_ci1(index(i)+1:index(i+1));
    Age_upper_ci(index(i)+1:index(i+1)) = Age_upper_ci1(index(i)+1:index(i+1));
    clear Intcal_A01 Age_const1 Age_lower_ci1 Age_upper_ci1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Checking for longer-term DCF trends to find obvious growth rate changes
    %%% if enough data points avalaibly and growth section long enough (>500a)
    %%% otherwise DCF variability might hide any decay and the fitting result 
    %%% is too strongly affected 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(depth(index(i):index(i+1)-1))>=2*NRCV & Age_const(index(i+1))-Age_const(index(i)+1)>key_breakpoint_age_diff
        breakpoint = SlopeBreakpoint(log_fm1(index(i)+1:index(i+1))',depth,gerade1, [index(i)+1 index(i+1)+1],NRCV);
        number_bp = length(breakpoint);
    else
        breakpoint = [ ];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% if no longer-term residual 14C trends could be observed (or too few 14C
    %%% measurements available), do nothing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(breakpoint)==1   
    
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% otherwise split the record into 'number_bp' growth sections and repeat 
    %%% the fitting procedure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        hold on
 
        for ii = 1 : number_bp
            plot([breakpoint(ii) breakpoint(ii)], [min(Age_const) max(Age_const)],'g')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            index_bp1(ii) = find(abs(depth-breakpoint(ii))==min(abs(depth-breakpoint(ii))) );
        end

        if depth(index(i)-1) <= depth(index_bp1(1))   
        %%% fitting the first part
            for jj = 1: index_bp1(1) - index(i)
                x3{jj} = x2{jj};
                y3{jj} = y2{jj};
            end
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR2_fit(depth, fm, [index(i) index_bp1(1)+1], Age_const(index(i)), depth(index(i)), x3, y3);

            Intcal_A0(index(i)+1:index_bp1(1)) = Intcal_A0_dummy(index(i)+1:index_bp1(1));
            Age_const(index(i)+1:index_bp1(1)) = Age_const_dummy(index(i)+1:index_bp1(1));
            Age_lower_ci(index(i)+1:index_bp1(1)) = Age_lower_ci_dummy(index(i)+1:index_bp1(1));
            Age_upper_ci(index(i)+1:index_bp1(1)) = Age_upper_ci_dummy(index(i)+1:index_bp1(1));
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

        %%% fitting the middle part
            for ii=2:number_bp
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(index_bp1(ii-1)) Age_upper_ci(index_bp1(ii-1)) Age_lower_ci(index_bp1(ii-1))], depth(index_bp1(ii-1)));

                Intcal_A0(index_bp1(ii-1)+1:index_bp1(ii)) = Intcal_A0_dummy(2:end);
                Age_const(index_bp1(ii-1)+1:index_bp1(ii)) = Age_const_dummy(2:end);
                Age_lower_ci(index_bp1(ii-1)+1:index_bp1(ii)) = Age_lower_ci_dummy(2:end);
                Age_upper_ci(index_bp1(ii-1)+1:index_bp1(ii)) = Age_upper_ci_dummy(2:end);
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            
        %%% fitting the oldest part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index_bp1(end) index(i+1)+1], [Age_const(index_bp1(end)) Age_upper_ci(index_bp1(end)) Age_lower_ci(index_bp1(end))], depth(index_bp1(end)));
    
            Intcal_A0(index_bp1(end)+1:index(i+1)) = Intcal_A0_dummy(2:end);
            Age_const(index_bp1(end)+1:index(i+1)) = Age_const_dummy(2:end);
            Age_lower_ci(index_bp1(end)+1:index(i+1)) = Age_lower_ci_dummy(2:end);
            Age_upper_ci(index_bp1(end)+1:index(i+1)) = Age_upper_ci_dummy(2:end);



        elseif depth(index(i)-1) >= depth(index_bp1(end)) 
        %%% fitting the oldest part
            for jj = 1: index(i+1) - index_bp1(end) +1
                x3{jj} = x2{size(x2,2)-((index(i+1))-index_bp1(end)+1)+jj};
                y3{jj} = y2{size(x2,2)-((index(i+1))-index_bp1(end)+1)+jj};
            end
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR2_fit(depth, fm, [index_bp1(end)-1 index(i+1)+1], Age_const(1), depth(1), x3, y3);

            Intcal_A0(index_bp1(end):index(i+1)) = Intcal_A0_dummy(Age_const_dummy > 0);
            Age_const(index_bp1(end):index(i+1)) = Age_const_dummy(Age_const_dummy > 0);
            Age_lower_ci(index_bp1(end):index(i+1)) = Age_lower_ci_dummy(Age_const_dummy > 0);
            Age_upper_ci(index_bp1(end):index(i+1)) = Age_upper_ci_dummy(Age_const_dummy > 0);
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

        %%% fitting the middle part
            for ii=number_bp:-1:2
                [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                    GR_fit(depth, fm, [index_bp1(ii-1) index_bp1(ii)+1], [Age_const(index_bp1(ii)) Age_upper_ci(index_bp1(ii)) Age_lower_ci(index_bp1(ii))], depth(index_bp1(ii)));

                Intcal_A0(index_bp1(ii-1):index_bp1(ii)-1) = Intcal_A0_dummy(1:end-1);
                Age_const(index_bp1(ii-1):index_bp1(ii)-1) = Age_const_dummy(1:end-1);
                Age_lower_ci(index_bp1(ii-1):index_bp1(ii)-1) = Age_lower_ci_dummy(1:end-1);
                Age_upper_ci(index_bp1(ii-1):index_bp1(ii)-1) = Age_upper_ci_dummy(1:end-1);
                clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy

            end
            
        %%% fitting the first part
            [log_fm1, gerade1, Intcal_A0_dummy, Age_const_dummy, Age_lower_ci_dummy, Age_upper_ci_dummy]  = ...
                GR_fit(depth, fm, [index(i)+1 index_bp1(1)+1], [Age_const(index_bp1(1)) Age_upper_ci(index_bp1(1)) Age_lower_ci(index_bp1(1))], depth(index_bp1(1)));

            Intcal_A0(index(i)+1:index_bp1(1)-1) = Intcal_A0_dummy(1:end-1);
            Age_const(index(i)+1:index_bp1(1)-1) = Age_const_dummy(1:end-1);
            Age_lower_ci(index(i)+1:index_bp1(1)-1) = Age_lower_ci_dummy(1:end-1);
            Age_upper_ci(index(i)+1:index_bp1(1)-1) = Age_upper_ci_dummy(1:end-1);
            clear Intcal_A0_dummy Age_const_dummy Age_lower_ci_dummy Age_upper_ci_dummy


  
            
        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
    
    



%%%Check for DCF changes in last stalagmite growth section

if isempty(i)==1, i = 1;end
%Just in case there is no hiatus, theen this line guarantees, that there will be some output in Fig3

lambda = 0.00012096;
stal_ainit = fm(index(i)+1:index(i+1)-1) ./ exp(-lambda*Age_const(index(i)+1:index(i+1)-1));
DCF = 1-(stal_ainit./Intcal_A0(index(i)+1:index(i+1)-1));% DCF of the top samples

avg_DCF_top = mean(DCF); %use this DCF to estimate the "starting point" of the second part (below the hiatus)
DCF_STD = std(DCF); %Assume constant DCF over the hiatus
fprintf('Average DCF from growth section %g: %g +/- %G %% \n', i,avg_DCF_top*100,DCF_STD*100)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate output files and most important figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if key<=4
    figure(2)
    hold on

    plot(depth(1:index(2)), Age_const(1:index(2)), 'b')
    plot(depth(1:index(2)), Age_lower_ci(1:index(2)), 'r')
    plot(depth(1:index(2)), Age_upper_ci(1:index(2)), 'r')
    for i = 1:length(index)-1
        plot(depth(index(i)+1:index(i+1)), Age_const(index(i)+1:index(i+1)), 'b')
        plot(depth(index(i)+1:index(i+1)), Age_lower_ci(index(i)+1:index(i+1)), 'r')
        plot(depth(index(i)+1:index(i+1)), Age_upper_ci(index(i)+1:index(i+1)), 'r')
    end
end

stal_ainit = fm ./ exp(-lambda*Age_const);
DCF = 1-(stal_ainit./Intcal_A0);% DCF of the top samples

if key == 5 % saving data only for real data sets

    fid = fopen([value{6}(1:end-4),'_results.txt'],'w');
    results(:,1) = depth;
    results(:,2) = fm;
    results(:,3) = fm_err;
    results(:,4) = Age_const;
    results(:,5) = Age_lower_ci;
    results(:,6) = Age_upper_ci;
    fprintf(fid,'depth\tfm\tfm_err\tage\tupper limit(1s)\tlower limit (1s)\r\n');
    fprintf(fid,'[mm]\t[fm]\t[fm]\t[a]\t[a]\t[a]\r\n');
    fprintf(fid,'%5.2f\t%5.4f\t%5.4f\t%5.0f\t%5.0f\t%5.0f\r\n', results');
    fclose(fid);
    clear results


    if strcmpi(proxyfile_d13C,'none') == 0        
        [depth1, d13C] = textread(proxyfile_d13C,...
            '%f %f','commentstyle','matlab');

        %%% interpolating radiocarbon time scale to proxy timescale and save
        %%% results: WITHOUT extrapolation

        for ii = 2:length(index)
            Age_const1 = interp1(depth(index(ii-1):index(ii)), Age_const(index(ii-1):index(ii)), depth1, 'linear');
            Age_lower_ci1 = interp1(depth(index(ii-1):index(ii)), Age_lower_ci(index(ii-1):index(ii)), depth1, 'linear');
            Age_upper_ci1 = interp1(depth(index(ii-1):index(ii)), Age_upper_ci(index(ii-1):index(ii)), depth1, 'linear');
        end

        fid = fopen([ value{7}(1:end-4),'_results.txt'],'w');
        results(:,1) = depth1;
        results(:,2) = d13C;
        results(:,3) = Age_const1;
        results(:,4) = Age_lower_ci1;
        results(:,5) = Age_upper_ci1;
        fprintf(fid,'depth\td13C\tage\tupper limit(1s)\tlower limit (1s)\r\n');
        fprintf(fid,'[mm]\t[permil]\t[a]\t[a]\t[a]\r\n');

        %%% determine significant digits for output files
        Age_const1_min = min(diff(Age_const1));
        if Age_const1_min < 1 & Age_const1_min >= 0.1
            fprintf(fid,'%5.2f\t%4.2f\t%5.1f\t%5.1f\t%5.1f\r\n', results');
        elseif Age_const1_min < 0.1 & Age_const1_min >= 0.01
            fprintf(fid,'%5.2f\t%4.2f\t%5.2f\t%5.2f\t%5.2f\r\n', results');
        elseif Age_const1_min < 0.01
            fprintf(fid,'%5.2f\t%4.2f\t%5.3f\t%5.3f\t%5.3f\r\n', results');
        else
            fprintf(fid,'%5.2f\t%4.2f\t%5.0f\t%5.0f\t%5.0f\r\n', results');
        end

        fclose(fid);
        clear results
    end

    if strcmpi(proxyfile_MgCa,'none') == 0        
        [depth2, MgCa] = textread(proxyfile_MgCa,...
            '%f %f','commentstyle','matlab');

        %%% interpolating radiocarbon time scale to proxy timescale and save
        %%% results: WITHOUT extrapolation

        for ii = 2:length(index)
            Age_const2 = interp1(depth(index(ii-1):index(ii)), Age_const(index(ii-1):index(ii)), depth2, 'linear');
            Age_lower_ci2 = interp1(depth(index(ii-1):index(ii)), Age_lower_ci(index(ii-1):index(ii)), depth2, 'linear');
            Age_upper_ci2 = interp1(depth(index(ii-1):index(ii)), Age_upper_ci(index(ii-1):index(ii)), depth2, 'linear');
        end

        fid = fopen([value{8}(1:end-4),'_results.txt'],'w');
        results(:,1) = depth2;
        results(:,2) = MgCa;
        results(:,3) = Age_const2;
        results(:,4) = Age_lower_ci2;
        results(:,5) = Age_upper_ci2;
        fprintf(fid,'depth\tMgCa\tage\tupper limit(1s)\tlower limit (1s)\r\n');
        fprintf(fid,'[mm]\t[ratio]\t[a]\t[a]\t[a]\r\n');
        %%% determine significant digits for output files
        Age_const2_min = min(diff(Age_const2));
        if Age_const2_min < 1 & Age_const2_min >= 0.1
            fprintf(fid,'%5.2f\t%5.4f\t%5.1f\t%5.1f\t%5.1f\r\n', results');
        elseif Age_const2_min < 0.1 & Age_const2_min >= 0.01
            fprintf(fid,'%5.2f\t%5.4f\t%5.2f\t%5.2f\t%5.2f\r\n', results');
        elseif Age_const2_min < 0.01
            fprintf(fid,'%5.2f\t%5.4f\t%5.3f\t%5.3f\t%5.3f\r\n', results');
        else
            fprintf(fid,'%5.2f\t%5.4f\t%5.0f\t%5.0f\t%5.0f\r\n', results');
        end

        fclose(fid);
        clear results
    end

end

figure(100)
subplot(2,5,[4,5])
box on
yyaxis right
plot(depth,DCF*100,'+-', 'color',[0.7 0.7 0.7])
ylabel('DCF [%]','Fontsize',14,'FontWeight','bold')
xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])


yyaxis left
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0],'YTick', [ ])

%subplot(2,5, [1,2,3,6,7,8])
subplot('Position', [0.1 0.11 0.45 0.815]);
hold on
box on

plot(depth(1:index(2)), Age_const(1:index(2)), 'color', [0.8 0.8 0.8],'LineWidth',2)
plot(depth(1:index(2)), Age_lower_ci(1:index(2)), ':', 'color', [0.8 0.8 0.8],'LineWidth',2)
plot(depth(1:index(2)), Age_upper_ci(1:index(2)), ':', 'color', [0.8 0.8 0.8],'LineWidth',2)
for i = 1:length(index)-1
    plot(depth(index(i)+1:index(i+1)), Age_const(index(i)+1:index(i+1)), 'color', [0.8 0.8 0.8],'LineWidth',2)
    plot(depth(index(i)+1:index(i+1)), Age_lower_ci(index(i)+1:index(i+1)), ':', 'color', [0.8 0.8 0.8],'LineWidth',2)
    plot(depth(index(i)+1:index(i+1)), Age_upper_ci(index(i)+1:index(i+1)), ':', 'color', [0.8 0.8 0.8],'LineWidth',2)
end

c(1:length(depth))=min(Age_const);
plot(depth,c,'^', 'color', [0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]) 
ylabel('Age [a BP]','Fontsize',14,'FontWeight','bold')
xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold')



