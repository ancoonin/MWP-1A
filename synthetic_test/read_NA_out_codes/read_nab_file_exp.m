%% read in convergence of expectation values from the nab output file

clear
close all

%-------------------------------------------------------------------------%                                                             
%%% This script reads the NA Bayes output file to check for convergence %%%
%%% of the expectation values for each free parameter                   %%%
%-------------------------------------------------------------------------%


fname='nab.out';


fid = fopen(fname);

il=0;

%find number of variables (set by nd in NA)

reading_on=1;

while reading_on==1
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>13 && isequal(tline(3:15),'Number of dim')
    C=strsplit(tline);
    nd=str2double(C{end}); %find number of variables
    reading_on=2; %i.e. move on to next task
    end
end


ranges=zeros(2,nd);
scales=zeros(nd,3);

% read in ranges of variables
while reading_on==2
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>=25 && strcmp(tline,'  Parameter space details :')

        %reading past next 4 lines to 1 after 'Parameter Ranges" line
        for i=1:5
        tline = fgetl(fid); 
        il=il+1;    
        end
        for i=1:nd
        tline = fgetl(fid); 
        il=il+1;    
        C=strsplit(tline);
        ind2remove=find(strlength(C)==0);
        C(ind2remove)=[];
        j=str2double(C(1)); %parameter #
        ranges(1,i)=str2double(C(2));
        ranges(2,i)=str2double(C(3));
        scales(i,:)=str2double(C(4:6));
        
            if j==nd
            reading_on=3;
            end
        end

    end
end

%find expected values in file
frewind(fid)
il = 0;
kp = 0;

while reading_on==3
tline = fgetl(fid); %reading current line
il=il+1;

if tline==-1
break
end

    if strlength(tline)>=9 && strcmp(tline(3:9),'Results')
    kp=kp+1;
    tline = fgetl(fid); %reading past blank line
    il=il+1;

    tline = fgetl(fid); %reading in line with number of random samples
    il=il+1;
    C=strsplit(tline);
    ind2remove=find(strlength(C)==0);
    C(ind2remove)=[];
    nmod(kp)=str2double(C(end));
    r=1;

        while r==1
        tline = fgetl(fid); %reading in expectation values
        il=il+1;
            if strlength(tline)>=30 && strcmp(tline(13:30),'Expectation values')
            tline = fgetl(fid); %reading past next two lines
            il=il+1;
    
            tline = fgetl(fid); 
            il=il+1;
                for i=1:nd
                tline = fgetl(fid); %reading in expectation values for parameter i
                il=il+1;
                newStr = split(tline,':');
                C=strsplit(newStr{2});
                ind2remove=find(strlength(C)==0);
                C(ind2remove)=[];
                data(1,kp,i)=str2double(C(1)); %mean expectation value
                data(2,kp,i)=str2double(C(2)); %standard error
                data(3,kp,i)=str2double(C(3)); %num error on that expectation value
                data(4,kp,i)=str2double(C(4)); %mean value for a uniform prior
                data(5,kp,i)=str2double(C(5)); %standard error for a uniform prior 


                r=2;
    
                end
            end
        end

    end
end

 %% now we will plot the convergence 
 m=ceil(sqrt(nd));
 n=m;
stitles={'EIS GMSL','WAIS GMSL','EIS onset','WAIS onset','EIS dur','WAIS dur'};
% plotting mean+/- numerical error
figure(2)
for p=1:nd
subplot(m,n,p)
mean=squeeze(data(1,:,p));
stnerr=squeeze(data(2,:,p));
numerr=squeeze(data(3,:,p));
meanplusne=mean+numerr;
meanminusne=mean-numerr;
meanpluserr=mean+stnerr;
meanminuserr=mean-stnerr;

plot(nmod,mean,'r-')
hold on
plot(nmod,meanplusne,'k--')
hold on
plot(nmod,meanminusne,'k--','HandleVisibility','off')

hold on
plot(nmod,meanpluserr,'b-.')
hold on
plot(nmod,meanminuserr,'b-.','HandleVisibility','off')
hold on
title(stitles{p})
xlabel('Ntot')
ylabel('Expectation Value')
    if p==nd
    legend({'mean','+/- numerical error','+/- standard error'},'Location','southoutside')
    end
end


%saving most numerically accurate expectation values for the mean value and
%standard error of each parameter
for p=1:nd
mean=squeeze(data(1,:,p));
stnerr=squeeze(data(2,:,p));
mean_final(p)=mean(end);
stnerr_final(p)=stnerr(end);
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',24)

end

 save('1Dstats.mat','mean_final',"stnerr_final")