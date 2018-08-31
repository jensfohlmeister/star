
function[log_fm1, gerade1, Intcal_A0, Age_const, Age_lower_ci, Age_upper_ci, lambda] = ...
    GR_fit(depth, fm, index, Anchorpoint_age, Anchorpoint_depth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Growth rate determination using the depth and radiocarbon content (fm)
% between two points (given by index(1) and index(2)) in the stalagmite and 
% the knowledge of an anchor point (anchorpoint_age with [mean min max] 
% and anchorpoint_depth) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read IntCal13, for later determination of atm. a14C values as new A0
[Intcal_age, ~, ~, Intcal_D14C, ~] = textread('intcal13.txt',...
'%f %f %f %f %f','commentstyle','matlab');

%%% Smooth the atmospheric 14C record to mimick the influence of vegetation
%%% and soil organic matter
windowsize = 23;
for i=(windowsize-1)/2+1:length(Intcal_D14C)-(windowsize-1)/2
    Intcal_D14C(i) = mean(Intcal_D14C(i-(windowsize-1)/2:i+(windowsize-1)/2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




iii = 1;
n=1; % stop next while conditions at n=21
Intcal_A0(:,1) = 1; %initial value, will be changed below
grate_1 = 1;  %growth_rate from earlier steps, now  1 mm/a

% The fit follows an equation of the type lnA/A0= -lambda*t
% Decay constant for 14C: lambda
lambda = 0.00012096;

%%% Age determination in section with anchor point %%%%%%%%%%%%%%%%%%%
while iii == 1

    %%%%%%%
    log_fm = log(fm(index(1):index(2)-1)./Intcal_A0);  
    %%%%%%%

    [fit,s]  = polyfit (depth(index(1):index(2)-1), log_fm,1);
    ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);


    % time is not known, but can be expressed as t=depth/growth rate
    % This results in lnA/A0 = -lambda*depth/growth rate
    % And lnA/A0*1/depth=slope
    slope = fit(:,1);
    grate = -lambda/slope;
    grate_error = grate .* sqrt((ste(1)./slope).^2);%Error on the growth rate (from error propagation)


    %%%%%%% Anchor point %%%%%
    Age_const = (depth(index(1):index(2)-1)-Anchorpoint_depth)/grate + Anchorpoint_age(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%


    for j=1:length(Age_const)
        idx1(j) =find(min(abs(Age_const(j)-Intcal_age))==abs(Age_const(j)-Intcal_age));
        Intcal_A0(j,1) = Intcal_D14C(idx1(j))/1000 + 1; % new A0 determined from IntCal13
    end
    
    if abs(grate-grate_1)<0.001
        iii = 0;
    else
        n=n+1;
        if n == 21
            iii=0;
            WRNGMSG = sprintf(['Attention! No good convergence of growth rate is achieved. \n' ...
                'Proceeding anyway with a close but not perfect fit. \n']);
            warning(WRNGMSG)
        end
        grate_1 = grate;
    end

end

%%% Age uncertainty (mainly due to error in slope fitting %%%%%%%%%%%%%%
for k = index(1): index(2)-1
    if depth(k)-Anchorpoint_depth <= 0
        Age_lower_ci(k-index(1)+1,1) = (depth(k)-Anchorpoint_depth)/abs(grate-grate_error) + Anchorpoint_age(3);
        Age_upper_ci(k-index(1)+1,1) = (depth(k)-Anchorpoint_depth)/(grate+grate_error) + Anchorpoint_age(2);
    else
        Age_lower_ci(k-index(1)+1,1) = (depth(k)-Anchorpoint_depth)/(grate+grate_error) + Anchorpoint_age(3);
        Age_upper_ci(k-index(1)+1,1) = (depth(k)-Anchorpoint_depth)/abs(grate-grate_error) + Anchorpoint_age(2);
    end
end
        
%%% Save results %%%%%%%%%%%%%%%%%%%
fprintf('mean growth rate for growth section between %g and %g mm: %f +- %f [mm/a] \n \n' ,...
    depth(index(1)), depth(index(2)-1), grate, grate_error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for age inversions (possible, when high degree of DCF variability
%%% compared to time window available for 14C decay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if grate<0
    ERRMSG = sprintf(['------------------------------------------------------------------- \n' ...
        'Age-depth modeling appears not appropriate in this section. Please, \n' ...
        'consider to set <key_breakpoint_age_diff> to at least 500years. The \n' ...
        'DCF variability seems to be relatively large for the short time \n' ...
        'interval covered by this stalagmite section.\n'...
        '------------------------------------------------------------------- \n']);
    error(ERRMSG)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% plot corrected ln(a14C) and the best fit line
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gerade1 = fit(1)*depth(index(1):index(2)-1)+fit(2);
log_fm1 = log(fm(index(1):index(2)-1)./Intcal_A0);                %corrected 14C values (atmospheric variations) 
% figure(1)
% plot(depth(index(1):index(2)-1),gerade1,'-r');
% plot(depth(index(1):index(2)-1),log_fm1,'+g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
