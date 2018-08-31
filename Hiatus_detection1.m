%%%%%%%%%Test for potential hiatuses in the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index,Anchorpoint_depth_index] = Hiatus_detection1(fm, depth, Anchorpoint_depth)

%%%%%%%%%We now use ln(fm) data to prevent
%%%%%%%%%the generation of too many false positive hiatuses.

log_fm = log(fm);

%Calculate difference between log(fm) values
depth_diff = diff(depth);
fm_diff = -diff(log_fm);
quotient_diff = fm_diff./depth_diff;

med = median(quotient_diff);
MAD = mad(quotient_diff); 


%Find fm differences that are more than 2.5MAD higher than the median (potential
%hiatus)
out= MAD*5;
clear index
index = find(quotient_diff > out);


%Plot the fm vs. depth relationship with numbers for the ages, to visualize
%and simplify choosing the samples for the users.
figure(6)
plot(depth, log_fm, 'o-k')
hold on
ylabel('ln(fM)','Fontsize',14,'FontWeight','bold')
xlabel('distance from top [mm]','Fontsize',14,'FontWeight','bold')
set(gca,'Fontsize',12,'LineWidth',2,'FontWeight','bold', 'YColor',[0 0 0])

numbers = [1: length(depth)]';
a = num2str(numbers);
b = cellstr(a);
dx = 0.01;
dy = 0.01;
text(depth+dx, log_fm+dy, b);
if isempty(index)
    str = sprintf('No growth stop detected. \nPlease confirm in command window.');
else
    str = sprintf('growth stop(s) probably at: %i \nPlease confirm in command window.', index);
end
legend(str);

%Ask whether the hiatuses should be incorporated
if isempty(index)==1
    clear index;
    prompt = 'No growth gaps detected: specify any points that should be considered a hiatus';
    x = input(prompt)
elseif isempty(index) ==0
    growth_gaps = index
    prompt = 'One or more growth gaps have been detected in the data (listed in growth_gaps). \n Please specify which points should be included in the model.';
     x= input(prompt)
end   


clear index
index = x;

%Corresponding points need to be identified in the original time series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% index marks the position of growth stops (plus start and end - for a
%%% smarter modeling later on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index1= index; %points chosen previously
clear index
index(1) = 1;  %add top index
index(2:length(index1)+1) = index1; %add points chosen by user
index(length(index1)+2) = length(fm);%+1; %add bottom index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In which part lies the anchor point?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = length(index):-1:2
    if Anchorpoint_depth <= depth(index(k))%-1)
        Anchorpoint_depth_index = k-1;
    end
end

