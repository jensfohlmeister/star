

function[m2,s2, mm1,ss1] = montecarlo(mittel,sigma)

[date,c14age, errC14age, d14C, errd14C] = textread(...
 'intcal13.txt','%f %f %f %f %f','commentstyle', 'matlab');

anzahl = 100;        % wieviel zufalls variablen normalverteilt werden sollen
                        % 10000 für gutes Bild
mittel = mittel*1000;%-59; %um auf 1950 zu sitzen
                                           %wichtig für intcal
sigma = sigma * 1000;



a14C=(1+d14C/1000)*100;
%tic


date_interp = 0:50000;
c14age_interp = interp1(date,c14age,0:50000);
d14C_interp = interp1(date,d14C,0:50000);
errd14C_interp = interp1(date,errd14C,0:50000);



for n=1:length(sigma)
    x = floor( mittel(n) + sigma(n) * randn(anzahl,1));
%    n
    for i = 1:anzahl

        if x(i) < 0 , x(i) = 0; end
        j=find(round(x(i)) == date_interp);
        y(i)=c14age_interp(j);
        yy(i)=d14C_interp(j);
        yy_err(i) = errd14C_interp(j);
        datedate(i)=date_interp(j);

    end

    s1=std(y);
    ss1(n)=sqrt(std(yy)^2 + (mean(yy_err))^2);
    m1=mean(y);
    mm1(n)=mean(yy);
    datedate1=mean(datedate);

    datedate2(n)=datedate1;
    mmm2=mm1(n)*exp(-datedate2(n)/8266);

    s2(n)=s1;
    m2(n)=m1;                       % m2 und mm2 gibt beides das 14C alter
    mm2(n)=log(mmm2/100)*(-8033);
end
%toc





% figure(1)
% clf
% hist(y,100)    
% figure(2)
% clf
% hist(x,100)


