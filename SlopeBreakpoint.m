
function[breakpoint] = SlopeBreakpoint(log_fm1,depth,gerade1,index,NRCV)

% checking the residual 14C variations for large scale signatures 
% especially if residuals can be well approximated by polynomials

significance3(2,1) = 1;
significance4(2,1) = 1;
significance5(2,1) = 1;
extrema_index3 = 1;
extrema_index4 = 1;
extrema_index5 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% checking for significance polynomial of order 2 to 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p2 = polyfit(depth(index(1):index(2)-1),log_fm1 - gerade1,2);
%%% Testing if polynomial fit of order 2 is good enough by correlation analysis
fitted_points = p2(3) + p2(2)*depth(index(1):index(2)-1) + p2(1) * depth(index(1):index(2)-1).^2;
[r, significance]= corrcoef(fitted_points, log_fm1 - gerade1);
% r(2,1);
% significance(2,1)

%%% to enable a later check if extreme values of polynomials are within the 
%%% length of the stal and not at the edge
not_at_edge = round(length(depth(index(1):index(2)-1))/10); %for finding extrema:disregarding the 10% at the edges


warning('OFF', 'all')
if length(log_fm1) >= 3*NRCV
    p3 = polyfit(depth(index(1):index(2)-1),log_fm1 - gerade1,3);
    %%% Testing if polynomial fit of order 3 is good enough by correlation analysis
    fitted_points3= p3(4) + p3(3)*depth(index(1):index(2)-1) + p3(2) * depth(index(1):index(2)-1).^2 ...
        + p3(1)*depth(index(1):index(2)-1).^3;
    [r3, significance3]= corrcoef(fitted_points3, log_fm1 - gerade1);
%     r3(2,1);
%     significance3(2,1)    
    diff_fitted_points3 = diff(fitted_points3);
    
    %%% to enable a later check if extreme values of polynomials are within the 
    %%% length of the stal and not at the edge
    not_at_edge = max(NRCV,round(length(depth(index(1):index(2)-1))/10)); %for finding extrema:disregarding the 10% at the edges

    extrema_index3 = 1+not_at_edge+find(diff_fitted_points3(not_at_edge+1:end-not_at_edge).*diff_fitted_points3(not_at_edge+2:end-not_at_edge+1)<0);
end
    
if length(log_fm1) >= 4*NRCV
    p4 = polyfit(depth(index(1):index(2)-1),log_fm1 - gerade1,4);
    %%% Testing if polynomial fit of order 4 is good enough by correlation analysis
    fitted_points4= p4(5) + p4(4)*depth(index(1):index(2)-1) + p4(3) * depth(index(1):index(2)-1).^2 ...
        + p4(2)*depth(index(1):index(2)-1).^3 + p4(1)*depth(index(1):index(2)-1).^4;
    [r4, significance4]= corrcoef(fitted_points4, log_fm1 - gerade1);
%     r4(2,1);
%     significance4(2,1)
    diff_fitted_points4 = diff(fitted_points4);
    
    %%% to enable a later check if extreme values of polynomials are within the 
    %%% length of the stal and not at the edge
    not_at_edge = max(NRCV,round(length(depth(index(1):index(2)-1))/10)); %for finding extrema:disregarding the 10% at the edges

    extrema_index4 = 1+not_at_edge+find(diff_fitted_points4(not_at_edge+1:end-not_at_edge).*diff_fitted_points4(not_at_edge+2:end-not_at_edge+1)<0);
end

if length(log_fm1) >= 5*NRCV
    p5 = polyfit(depth(index(1):index(2)-1),log_fm1 - gerade1,5);
    %%% Testing if polynomial fit of order 5 is good enough by correlation analysis
    fitted_points5= p5(6) + p5(5)*depth(index(1):index(2)-1) + p5(4) * depth(index(1):index(2)-1).^2 ...
        + p5(3)*depth(index(1):index(2)-1).^3 + p5(2)*depth(index(1):index(2)-1).^4 + p5(1)*depth(index(1):index(2)-1).^5;
    [r5, significance5]= corrcoef(fitted_points5, log_fm1 - gerade1);
%     r5(2,1);
%     significance5(2,1)
    diff_fitted_points5 = diff(fitted_points5);

    %%% to enable a later check if extreme values of polynomials are within the 
    %%% length of the stal and not at the edge
    not_at_edge = max(NRCV,round(length(depth(index(1):index(2)-1))/10)); %for finding extrema:disregarding the 10% at the edges

    extrema_index5 = 1+not_at_edge+find(diff_fitted_points5(not_at_edge+1:end-not_at_edge).*diff_fitted_points5(not_at_edge+2:end-not_at_edge+1)<0);
end
warning('ON', 'all')

polyorder=2;
if significance3(2,1)<0.05, polyorder=3; end
if significance4(2,1)<0.05, polyorder=4; end
if significance5(2,1)<0.05, polyorder=5; end

if polyorder==5
    if length(extrema_index5) < 4 | significance5(2,1)/significance4(2,1)>0.1, polyorder=4; end
end
if polyorder==4
    if length(extrema_index4) < 3 | significance4(2,1)/significance3(2,1)>0.1, polyorder=3; end
end
if polyorder==3
    if length(extrema_index3) < 2 | significance3(2,1)/significance(2,1)>0.1, polyorder=2; end
end
    
if polyorder==5, significance(2,1)=significance5(2,1);fitted_points=fitted_points5; end
if polyorder==4, significance(2,1)=significance4(2,1);fitted_points=fitted_points4; end
if polyorder==3, significance(2,1)=significance3(2,1);fitted_points=fitted_points3; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% two options to determine the breakpoint in growthrate
%%% 1. depth with DCF extreme value of fit (first derivative of  polynomial of order
%%%     two [y=a+bx+cx^2 -> y'=b+2cx] with y' = 0 --> solve for x)
%%% 2. depth with DCF extreme value of fit after smoothing (to get rid of
%%%     short-term DCF variations

if significance(2,1)>=0.05
    fprintf('After testing the residuals of 14C, it appears unlikely that a major shift in growthrate occured.\n \n')
    breakpoint =[ ];
else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% smoothing by a running mean of  5 points 
    %%% smoothing is neccessary to get rid of short-term DCF variations; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    smoothed_residual(1) = sum(log_fm1(1:3) - gerade1(1:3))/3;
    smoothed_residual(2) = sum(log_fm1(1:4) - gerade1(1:4))/4;
    for i = 3:length(log_fm1)-2
        smoothed_residual(i) = sum(log_fm1(i-2:i+2) - gerade1(i-2:i+2))/5;
    end
    smoothed_residual(i+1) = sum(log_fm1(i-1:i+2) - gerade1(i-1:i+2))/4;    
    smoothed_residual(i+2) = sum(log_fm1(i:i+2) - gerade1(i:i+2))/3;

    if polyorder == 2
        if p2(1)<0
            breakpoint = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(not_at_edge:end-not_at_edge+1))));
        else
            breakpoint = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(not_at_edge:end-not_at_edge+1))));
        end
    elseif polyorder == 3
        %%% search for extreme point within +/- 25% of all measurements around
        %%% the extreme points of the polynom
        low_lim(1) = extrema_index3(1) - round(length(log_fm1)/4);
        upp_lim(1) = extrema_index3(1) + round(length(log_fm1)/4);
        low_lim(2) = extrema_index3(2) - round(length(log_fm1)/4);
        upp_lim(2) = extrema_index3(2) + round(length(log_fm1)/4);
        for i = 1:2
            if low_lim(i) < not_at_edge, low_lim(i) = not_at_edge; end
            if upp_lim(i) > length(log_fm1)-not_at_edge+1, upp_lim(i) = length(log_fm1)-not_at_edge+1; end
        end
        
        if fitted_points3(extrema_index3(1)) > fitted_points3(extrema_index3(2))
            breakpoint(1) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(1):upp_lim(1)))));
            breakpoint(2) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(2):upp_lim(2)))));
        else
            breakpoint(1) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(1):upp_lim(1)))));
            breakpoint(2) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(2):upp_lim(2)))));
        end

        breakpoint = sort(breakpoint);
        
    elseif polyorder == 4
        %%% search for extreme point within +/- 25% of all measurements around
        %%% the extreme points of the polynom
        low_lim(1) = extrema_index4(1) - round(length(log_fm1)/5);
        upp_lim(1) = extrema_index4(1) + round(length(log_fm1)/5);
        low_lim(2) = extrema_index4(2) - round(length(log_fm1)/5);
        upp_lim(2) = extrema_index4(2) + round(length(log_fm1)/5);
        low_lim(3) = extrema_index4(3) - round(length(log_fm1)/5);
        upp_lim(3) = extrema_index4(3) + round(length(log_fm1)/5);

        for i = 1:3
            if low_lim(i) < not_at_edge, low_lim(i) = not_at_edge; end
            if upp_lim(i) > length(log_fm1)-not_at_edge+1, upp_lim(i) = length(log_fm1)-not_at_edge+1; end
        end
        
        if fitted_points4(extrema_index4(1)) > fitted_points4(extrema_index4(2))
            if fitted_points4(extrema_index4(2)) < fitted_points4(extrema_index4(3))
                breakpoint(1) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(1):upp_lim(1)))));
                breakpoint(2) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(2):upp_lim(2)))));
                breakpoint(3) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(3):upp_lim(3)))));
            else
                breakpoint(1) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(1):upp_lim(1)))));
                breakpoint(2) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(2):upp_lim(2)))));
                breakpoint(3) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(3):upp_lim(3)))));              
            end
        else
            if fitted_points4(extrema_index4(2)) > fitted_points4(extrema_index4(3))
                breakpoint(1) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(1):upp_lim(1)))));
                breakpoint(2) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(2):upp_lim(2)))));
                breakpoint(3) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(3):upp_lim(3)))));
            else
                breakpoint(1) = depth(index(1)-1+find(smoothed_residual==min(smoothed_residual(low_lim(1):upp_lim(1)))));
                breakpoint(2) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(2):upp_lim(2)))));
                breakpoint(3) = depth(index(1)-1+find(smoothed_residual==max(smoothed_residual(low_lim(3):upp_lim(3)))));

            end
        end

        breakpoint = sort(breakpoint);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% if no age depth modelling is possible (GR negative(!) due to too 
    %%% short period of time; too large DCF variations) with choosen 
    %%% breakpoints, shift breakpoints slightly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1 : length(breakpoint)
        index_bp1(ii) = find(abs(depth-breakpoint(ii))==min(abs(depth-breakpoint(ii))) ); 
    end
    
    [fit]  = polyfit (depth(index(1):index_bp1(1)), log_fm1(1:index_bp1(1)-index(1)+1),1);
    while fit(:,1)>0
        index_bp1(1) = index_bp1(1)+1;
        breakpoint(1) = depth(index_bp1(1));
        [fit]  = polyfit (depth(index(1):index_bp1(1)), log_fm1(1:index_bp1(1))-index(1)+1,1);
    end
    
    for ii = 2:length(breakpoint)
        [fit]  = polyfit (depth(index_bp1(ii-1):index_bp1(ii)), log_fm1(index_bp1(ii-1)-index(1)+1:index_bp1(ii)-index(1)+1),1);
    while fit(:,1)>0
        index_bp1(ii) = index_bp1(ii)+1;
        breakpoint(ii) = depth(index_bp1(ii));
        [fit]  = polyfit (depth(index_bp1(ii-1):index_bp1(ii)), log_fm1(index_bp1(ii-1)-index(1)+1:index_bp1(ii)-index(1)+1),1);
    end

    end
    
    [fit]  = polyfit (depth(index_bp1(end):end), log_fm1(index_bp1(end)-index(1)+1:end),1);
    while fit(:,1)>0
        index_bp1(end) = index_bp1(end)-1;
        breakpoint(end) = depth(index_bp1(end));
        [fit]  = polyfit (depth(index_bp1(end):end), log_fm1(index_bp1(end)-index(1)+1:end),1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    fprintf('It seems there is one major change in growth rate around %3.1f mm depth \n \n', breakpoint)
    fprintf('Recalculating growth rate with the assumption of phases with different growth rates. \n \n')

    figure(20)
    hold on
    box on
    ylabel('DCF offset [A.U.]','Fontsize',14,'FontWeight','bold')
    xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
    set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])

    plot(depth(index(1):index(2)-1),log_fm1 - gerade1, '*-k')
    plot(depth(index(1):index(2)-1),smoothed_residual, '*r')
    plot(depth(index(1):index(2)-1), fitted_points, '-b')
    plot([breakpoint(1) breakpoint(1)], [max(smoothed_residual) min(smoothed_residual)], '-g')
    if polyorder == 3
        plot([breakpoint(2) breakpoint(2)], [max(smoothed_residual) min(smoothed_residual)], '-g')
    end
    if polyorder == 4
        plot([breakpoint(2) breakpoint(2)], [max(smoothed_residual) min(smoothed_residual)], '-g')
        plot([breakpoint(3) breakpoint(3)], [max(smoothed_residual) min(smoothed_residual)], '-g')
    end

    legend('residual 14C','smoothed residual 14C', 'fitted polynom', sprintf( '%s\n%s', 'estimated best choice', 'of growth rate change' )' ,...
        'Location','best')

end















