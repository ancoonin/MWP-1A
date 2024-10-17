%% plot marginal posterior PDFs from nab.out

%------------------------------------------------------------------------%                                                             
%%% This script reads nab.out from the NA Bayes run and plots the      %%%
%%% marginal posterior probability distribution functions              %%%
%%% for the free parameters in the synthetic inversion test            %%%
%%% e.g., contribution to GMSL, onset time of melting, duration of     %%%
%%% melting for both th Eurasian and West Antarctic ice sheets         %%%
%------------------------------------------------------------------------%


clear
close all

icesheets={'Eurasian Ice Sheet', 'West Antarctic Ice Sheet'};
fname='nab.out';
fid = fopen(fname);

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

while reading_on==2
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>12 && isequal(tline(3:14),'Calculate 1D')
        if isequal(tline(end-2),'y')
            for i=1:3
            tline = fgetl(fid); %reading current line
            il=il+1;
            end
            tline = fgetl(fid); %reading current line
            il=il+1;
            binstr=tline(48:end);
            % number of bins per axis for 1D marginals
            ndis=str2double(binstr);
            np=nd; %number of plots
            reading_on=3; %moving on to next task
        end
    end
end

%find final results in file 
while reading_on==3
tline = fgetl(fid); %reading current line
il=il+1;
    if strcmp(tline,'  Results of Monte Carlo integration using NA random walk')
    ipos = il;
    reading_on=4;
    end
end

% find position in file of final 1D marginals
kp = 0;
frewind(fid)
for i=1:ipos
tline = fgetl(fid);
end

while reading_on==4
tline = fgetl(fid);
il=il+1;

if tline==-1
break
end
    if strlength(tline)>=10
    strcheck=tline(1:10);
        if strcmp(strcheck,'  Marginal')
        %read in 1-D marginals
            for i=1:np
            tline = fgetl(fid);
            il=il+1;
                for j=1:ndis 
                tline = fgetl(fid);
                il=il+1;
                C=strsplit(tline);
                ind2remove=find(strlength(C)==0);
                C(ind2remove)=[];
                data(j,1,i)=str2double(C(1)); %for each parameter i there are j values
                data(j,2,i)=str2double(C(2)); %each associated with a probability value

                end
                tline = fgetl(fid); 
                il=il+1;
                tline = fgetl(fid); 
                il=il+1;

                if i==nd
                kp = kp + 1;
                end
            end
        end
    end
end



%check if the pdfs are already normalized between 0 and 1:

for p=1:nd
A = trapz(data(:,1,p),data(:,2,p));

normfreq(:,:,p)=[squeeze(data(:,1,p)) squeeze(data(:,2,p))/A]; %probability distribution is the density distribution normalized by the integral of the density distribution (so the probability distribution integrates to 1)

Aprime(p)=trapz(normfreq(:,1),normfreq(:,2));

%normalized posterior PDFs (each integrates to 1)

end


%what is the maximum yval for the normalized probabilities?
ymaxindiv=max(normfreq(:,2,1:nd/3));
ymaxplot(1)=max(ymaxindiv);

%plot 1-d Marginals

m=nd/2; % # of rows in subplot
n=2; % # of columns in subplot

figure(1)
for p=1:np/3
subplot(m,n,p)
area(normfreq(:,1,p),normfreq(:,2,p))
yvals=linspace(0,ymaxplot(1),10);
stitle=icesheets{p};



xlim([0.1 10])
ylim([0 ymaxplot(1)])
yticks([0 ymaxplot(1)])
yticklabels([0 1])
xticks([1 3 5 7 9])

end


ymaxindiv=max(normfreq(:,2,nd/3+1:nd/3*2));
ymaxplot(2)=max(ymaxindiv);


for p=np/3+1:np/3*2
subplot(m,n,p)
area(-1*normfreq(:,1,p)/1e3,normfreq(:,2,p))
yvals=linspace(0,ymaxplot(2),10);
set(gca,'XDir','reverse')
ylim([0 ymaxplot(2)])
yticks([0 ymaxplot(2)])
yticklabels([0 1])
xlim([14 14.675])
xticks([14 14.2 14.4 14.6])

end

%
hold off
% 
ymaxindiv=max(normfreq(:,2,2*nd/3+1:nd));
ymaxplot(3)=max(ymaxindiv);

%duration 
for p=2*np/3+1:np
subplot(m,n,p)
area(normfreq(:,1,p),normfreq(:,2,p))
yvals=linspace(0,ymaxplot(3),10);
ylim([0 ymaxplot(3)])
yticks([0 ymaxplot(3)])
yticklabels([0 1])
xticks([0 250 500])

end

%find modal values for each parameter
for p=1:np
modeid=find(normfreq(:,2,p)==max(normfreq(:,2,p)));
modeval(p)=normfreq(modeid,1,p); 
end


set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.04 0.1])
set(findall(gcf,'-property','Box'),'Box','on')


