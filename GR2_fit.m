function[log_fm2, gerade2, Intcal_A0, Age_const, Age_lower_ci, Age_upper_ci, lambda] = ...
    GR2_fit(depth, fm, index, Anchorpoint_age, Anchorpoint_depth, x2, y2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Growth rate determination using the depth and radiocarbon content (fm)
% between two points (given by index(1) and index(2)) in the stalagmite, 
% where tha anchor point is replaced by the known mean DCF and calibration
% of the radiocarbon values (x2, y2).
% The knowledge of a precise anchor point (anchorpoint_age with [mean] 
% and anchorpoint_depth) is not important, but a necessary for the
% computation process. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% read IntCal13, for determination of atm. a14C values at defined
    %%% ages and calculation of stalagmite 14C value
    [Intcal_age, ~, ~, Intcal_D14C, ~] = textread('intcal13.txt',...
    '%f %f %f %f %f','commentstyle','matlab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Smooth the atmospheric 14C record to mimick the influence of vegetation
%%% and soil organic matter
windowsize = 23;
for i=(windowsize-1)/2+1:length(Intcal_D14C)-(windowsize-1)/2
    Intcal_D14C(i) = mean(Intcal_D14C(i-(windowsize-1)/2:i+(windowsize-1)/2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Intcal_A0(index(1)+1:index(2)-1,1) = 1; %first guess, will be iteratively determined
grate_1 = 1;  %first guess of growth rate Wachstumsrate, will be iteratively determined
    
iij = 1;
n=0;
while iij == 1

    %%%%%%% log of the atmospheric 14C corrected values
    log_fm(index(1)+1:index(2)-1) = log(fm(index(1)+1:index(2)-1)./Intcal_A0(index(1)+1:index(2)-1));  
    log_fm = reshape(log_fm, [], 1);
    %%%%%%%
        
    %%% fitting the data %%%%%%%%%%%%%%%%%%%
    [fit2,s] = polyfit (depth(index(1)+1:index(2)-1), log_fm(index(1)+1:index(2)-1),1);
    ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);

    % The fit follows an equation of the type lnA/A0= -lambda*t
    % Decay constant for 14C: lambda
    lambda = 0.00012096;


    % time is not known, but can be expressed as t=depth/growth rate
    % This results in ln(A/A0) = -lambda*depth/growth rate
    % therefore the slope of the ln(A/A0)in dependence of depth is equal to -lamda/growthrate
    slope = fit2(:,1);
    grate = -lambda./slope;
    grate_error = grate .* sqrt((ste(1)./slope).^2);%Error on the growth rate (from error propagation)
 
    %%%%%%%
    % Estimate rough position of constant growth line, which  is needed to find
    % the youngest position of constanst growth line
    % for that purpose use the oldest age of the young section
    Age_const(index(1)+1:index(2)-1) = (depth(index(1)+1:index(2)-1)-Anchorpoint_depth)/grate + Anchorpoint_age(1);
    %%%%%%%

    %%% find the youngest position of all calibrated ages with respect to the
    %%% slope and put the growth rate through that point 
    for ii = index(1)+1:index(2)-1
        minimum(ii) = x2{ii-index(1)}(end)-Age_const(ii);
    end
    imin=find(minimum==min(minimum(index(1)+1:index(2)-1)));
    Age_const(index(1)+1:index(2)-1) = ...
        (depth(index(1)+1:index(2)-1)-depth(imin))/grate + x2{imin-index(1)}(end); %youngest line

    %%% shift the slope through all 14C prob densities and calculate the
    %%% weighted maximum for best placement
    count=1;
    sum_prob(count)=0.00001;

    while sum_prob(count) > 0 
    
        count = count+1;
        Age_const_shift(index(1)+1:index(2)-1) = Age_const(index(1)+1:index(2)-1)+count; % shift the line
        sum_prob(count)=0;

        for ii = index(1)+1:index(2)-1
            %%%Maybe the following 'if condition' could be speeded up somehow ...
            if min(abs(Age_const_shift(ii)-x2{ii-index(1)}))<=0.5
                sum_prob(count) = sum_prob(count)+...
                    y2{ii-index(1)}(find(abs(Age_const_shift(ii)-x2{ii-index(1)}(:))==min(abs(Age_const_shift(ii)-x2{ii-index(1)})))) - ...
                    depth(ii);  % subtract depth; was previously just added for 
                                % illustration purposes (for Fig. 10)
            end
        end

    end

    %%% find position of highest cumulative probability which is the best
    %%% estimate for the timing of the stalgmite growth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% the better way seems to be using the weighted mean of the probabillity
    %%% distribution (sum_prob), which is previously computed by the summed
    %%% probabilities at the intersection point of the growth line and the 
    %%% individual calibrated ages
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sum_prob = sum_prob/sum(sum_prob);
    count_vector = 1:length(sum_prob);
    kkk = sum(sum_prob.*count_vector);   %%% mean shift with respect to youngest 
                                 %%% position of growth slope
    %%% symetrical error
    kkk_err = sqrt(sum((sum_prob.*(count_vector-kkk)).*(count_vector-kkk)));
                                 %%% 1sigma error of position of slope
       
    Age_const(index(1)+1:index(2)-1) = Age_const(index(1)+1:index(2)-1) + kkk;

    for j=index(1)+1:index(2)-1
        idx2(j-index(1)) =find(min(abs(Age_const(j)-Intcal_age))==abs(Age_const(j)-Intcal_age));
        Intcal_A0(j) = Intcal_D14C(idx2(j-index(1)))/1000 + 1; % new A0 determined from IntCal13
    end                             

    if abs(grate-grate_1)<0.001
        iij = 0;
    else
        n=n+1;
        if n == 21
            iij=0;
            WRNGMSG = sprintf(['Attention! No good convergence of growth rate is achieved. \n' ...
                'Proceeding anyway with a close but not perfect fit. \n']);
            warning(WRNGMSG)
        end
        grate_1 = grate;
    end
        
%     figure(12)
%     plot(1:length(sum_prob),sum_prob)
%     hold on 
%     plot([kkk kkk], [0 max(sum_prob)],'-k') % Mean of probability distribution
%                                         % and used as anchorpoint
end

%%% Age uncertainty (mainly due to error in slope fitting %%%%%%%%%%%%%%
for k = index(1)+1: index(2)-1
    if depth(k)-Anchorpoint_depth <= 0
        Age_lower_ci(k-index(1)+1,1) = (depth(k)-depth(index(2)-1))/abs(grate-grate_error) + Age_const(index(2)-1);
        Age_upper_ci(k-index(1)+1,1) = (depth(k)-depth(index(2)-1))/(grate+grate_error) + Age_const(index(2)-1);
    else
        Age_lower_ci(k-index(1)+1,1) = (depth(k)-depth(index(1)+1))/(grate+grate_error) + Age_const(index(1)+1);
        Age_upper_ci(k-index(1)+1,1) = (depth(k)-depth(index(1)+1))/abs(grate-grate_error) + Age_const(index(1)+1);
    end
end
%%% plus age error due to uncertainty in DCF
Age_lower_ci(index(1)+1:index(2)-1) = Age_lower_ci(2:end) - kkk_err;
Age_upper_ci(index(1)+1:index(2)-1) = Age_upper_ci(2:end) + kkk_err;
Age_lower_ci(1:index(1)) = 0;
Age_upper_ci(1:index(1)) = 0;


%%% Save results %%%%%%%%%%%%%%%%%%%
grate_save(i) = grate;
grate_error_save(i) = grate_error;

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



% figure(13) % plotting the best age-depth model with 14C probabilities
% plot (Age_const(index(1)+1:index(2)-1), depth(index(1)+1:index(2)-1), 'r');
% plot (Age_lower_ci(index(1)+1:index(2)-1), depth(index(1)+1:index(2)-1), 'g');
% plot (Age_upper_ci(index(1)+1:index(2)-1), depth(index(1)+1:index(2)-1), 'g');
% 
% figure(1)
% % Age_lower = (depth(index:end)-depth(14))/(grate2+grate_error2) + Age_const2(1) - 2*kkk_err;
% % Age_upper = (depth(index:end)-depth(14))/(grate2-grate_error2) + Age_const2(1) + 2*kkk_err;
gerade2 = fit2(1)*depth(index(1)+1:index(2)-1)+fit2(2);
log_fm2(index(1)+1:index(2)-1) = log(fm(index(1)+1:index(2)-1)./Intcal_A0(index(1)+1:index(2)-1)); %14C increase in the past
% 
% plot(depth(index(1)+1:index(2)-1),gerade2,'-r');
% plot(depth(index(1)+1:index(2)-1),log_fm2(index(1)+1:index(2)-1),'+g');
% 
clear Age_corr Age_corr_err

