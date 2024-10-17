%% plot correlation/covariance matrix from nab output file

clear
close all
%-------------------------------------------------------------------------%                                                             
%%% This script reads in the posterior covariances (output of NA Bayes) %%%
%%% between the free parameters in the inversion.                       %%%
%-------------------------------------------------------------------------%


fname='nab.out';

fid = fopen(fname);
numIce=2;

reading_on=1;
il=0; %index of the line we are reading in the file
while reading_on==1
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>13 && isequal(tline(3:15),'Number of dim')
    C=strsplit(tline);
    nd=str2double(C{end}); %find number of variables
    reading_on=2; %i.e. move on to next task
    end
end

%find final results in file 
while reading_on==2
tline = fgetl(fid); %reading current line
il=il+1;
    if strcmp(tline,'  Results of Monte Carlo integration using NA random walk')
    ipos = il;
    reading_on=3;
    end
end

% find position in file of results
kp = 0;
frewind(fid)
for i=1:ipos
tline = fgetl(fid);
end

% read in final covariance matrix
while reading_on==3
tline = fgetl(fid);
il=il+1;

if tline==-1
break
end

    if strlength(tline)>=12 && strcmp(tline(3:12),'Covariance')
    kp = kp + 1;
    tline = fgetl(fid); %reading past blank line
    il=il+1;

        for i=1:nd
        tline = fgetl(fid); %reading past descriptive line
        il=il+1;
        %read in row of covariance matrix:
        %                   NAB.out is formatted so that there are 5 values in each line
        %                   this means we need to read in the 3 rows that give the values of 
        %                   the covariance matrix for the row corresponding to parameter i
            for k=1:ceil(nd/5)
            tline = fgetl(fid);
            il=il+1;
            C=strsplit(tline);
            ind2remove=find(strlength(C)==0);
            C(ind2remove)=[];
                if k==1
                covlist=str2double(C);
                else
                covlist=[covlist str2double(C)];
                end

            end
        
            cov(i,1:nd)=covlist;

        end 
    
        %reading past the next 3 lines to get to the numerical error in the cov matrix        
        for i=1:3
        tline = fgetl(fid);
        il=il+1;
        end
    
        %reading in cov matrix numerical error 
        for i=1:nd
        tline = fgetl(fid); %header
        il=il+1;
        %read in row of covariance matrix err


            for k=1:ceil(nd/5)
            tline = fgetl(fid);
            il=il+1;
            C=strsplit(tline);
            ind2remove=find(strlength(C)==0);
            C(ind2remove)=[];

                if k==1
                coverrlist=str2double(C);
                else
                coverrlist=[coverrlist str2double(C)];
                end

            end
     
          coverr(i,:)=coverrlist;

        end 

    end

end


%computing standard deviation as square root of the variance (diagonal entries of the covariance matrix)

for i=1:nd
sigma(i)=sqrt(cov(i,i));
end

save('standard_dev.mat','sigma');

icesheetDOF={'EIS GMSL','WAIS GMSL','EIS onset','WAIS onset','EIS duration','WAIS duration'};

c = redblue_cmap(11);

%% correlation matrix 
% For the off-diagonal elements we can plot the correlation matrix. 
% (since it is difficult to display the covariance matrix directly, given
% the parameters have different dimensions


for i=1:nd
    for j=1:nd
        gamma(i,j)=cov(i,j)/sqrt(cov(i,i)*cov(j,j));
    end
end


ii = ones(size(gamma));
idx = tril(ii);
gamma(~idx)=nan;

figure()
heatmap(gamma)
colormap(c)
colorbar
title('Correlation Matrix')
clim([-1 1])
ax = gca;


warning('off','MATLAB:structOnObject')
axp = struct(ax);     
axp.Axes.XAxisLocation = 'top';


set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.04 0.1])
set(findall(gcf,'-property','TickDir'),'TickDir','out'); 
set(findall(gcf,'-property','Box'),'Box','on')



priorcov=zeros(5); 


%take prior to be the difference between the min and max range ^2 

%    EIS      WAIS
%    0         0      
%    7.2       9.1  


 
prior(1,1)= 7.2;
prior(2,2)= 9.1;
prior(3,3)= 650; %onset
prior(4,4)= 650; %onset
prior(5,5)= 500; %onset
prior(6,6)= 500; %onset

invprior=1./(prior.^2).*eye(nd);
invprior(isinf(invprior))=0;

Resolution=eye(nd)-invprior.*cov;


% non dimensional resolution matrix

sig=sqrt(diag(prior)); %sigma i is a representative scale length for parameter i.
for i=1:nd
   
    for j=1:nd
        Rprime(i,j)=Resolution(i,j)*sig(j)/sig(i);
    end
end

figure()
heatmap(Rprime)
colormap(hot)
colorbar
title('Non-dimensional Resolution Matrix')
ax = gca;
clim([0 1])

warning('off','MATLAB:structOnObject')
axp = struct(ax);     
axp.Axes.XAxisLocation = 'top';


set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.04 0.1])
set(findall(gcf,'-property','TickDir'),'TickDir','out'); 
set(findall(gcf,'-property','Box'),'Box','on')






