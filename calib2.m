

function[x2, y2] = calib2(x, x_err,depth)
                        % ( 14C age, 14C age error, cm)

%for calibration purposes of all DCF corrected 14C values and a nice plot                        


[date1, C14Age1, C14Age_err1, ~, ~] = textread(...
    'intcal13.txt', '%f %f %f %f %f','commentstyle', 'matlab'); 
%%%annual interpolation
date = (50000:-1:0)';
%keyboard
C14Age = interp1(date1,C14Age1,date);
C14Age_err = interp1(date1,C14Age_err1,date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i = date(end-1):date(1)
%    C14Age(i) = -8033*log(aC14atm(i)/100*exp(-date(i)/8267));
%    C14Age_err(i) = (8033)/aC14atm(i)/100/exp(-(date(i)-56)/8267)*100*exp(-(date(i)-56)/8267)*aC14atm_err(i);
%end


for j = 1:length(x)
    %%% calibrate
    sigma1 = sqrt(C14Age_err.^2+x_err(j)^2);
    value = normpdf(x(j), C14Age, sigma1);

    value = value/sum(value);

    %%% removing unneccessary data
%     for i = 1:length(value)
% 
%         if value(i)/max(value)<0.001, value(i) = 0; end
%     end
    
    istart = find(value/max(value) >=0.001,1,'first');
    iend = find(value/max(value)>=0.001,1,'last');
    value(1:istart-1) = 0;
    value(iend+1:end) = 0;

    y2{j} = value(value.*date(1:end) ~= 0)';
    x2{j} = date(value.*date(1:end) ~= 0); 
end


length_of_stal = depth(end)+depth(end)/10;
max_height = length(depth)/length_of_stal; % number of samples by length of stal


% figure(10)
% hold on
% set(gca,'YLim', [0 length_of_stal+1/max_height], 'XDir','reverse',...
%     'Fontsize',12,'LineWidth',2,'FontWeight','bold')
% ylabel(['distance from top [cm]'],'Fontsize',14,'FontWeight','bold') 
% xlabel(['cal age [a b1950]'],'Fontsize',14,'FontWeight','bold') 

for j = 1:length(x)
    %%% scale values to sample depth
    scale(j) = max(y2{j});
    y2{j} = y2{j}/scale(j)/max_height;

    %%% shift values to according sample depth
    y2{j} = y2{j} + depth(j);
    
%     plot(x2{j}, y2{j})
%     plot([x2{j}(1) x2{j}(end)],[depth(j) depth(j)]) 

end







