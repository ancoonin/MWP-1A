%% Calculate 2D marginal posterior PDFs from nab.out 

%------------------------------------------------------------------------%                                                            
%%% This script reads nab.out from the NA Bayes run and plots the      %%%
%%% contour map of the joint posterior probability distribution        %%%
%%% for parameter choices set in nab.in (see marg2dkey array)          %%%
%%%                                                                    %%%
%%% The number and type of plots can be adjusted in nab.in             %%%
%%% see NA Source code) for more info on nab.in                        %%% 
%------------------------------------------------------------------------%

clear
close all
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

% find number of 2-D marginal plots to make

while reading_on==2
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>12 && isequal(tline(3:14),'Number of 2D')
    numplotstr=tline(48:end);
    np=str2double(numplotstr);
    reading_on=3; %move on to new task
    end
end

% find bin discretization 

while reading_on==3
tline = fgetl(fid); %reading current line
il=il+1;
    if strlength(tline)>=32 && isequal(tline(31:32),'2D')
    binstr=tline(48:end);
    % number of bins per axis for 2D marginals
    ndis=str2double(binstr);   
    reading_on=4; %move on to next task
    end
end

%find final results in file
while reading_on==4
tline = fgetl(fid); %reading current line
il=il+1;
    if strcmp(tline,'  Results of Monte Carlo integration using NA random walk')
    ipos = il;
    reading_on=5;
    end
end

% final 2D marginals
kp = 0;
frewind(fid)
for i=1:ipos
tline = fgetl(fid);
end

while reading_on==5

tline = fgetl(fid);
il=il+1;

if tline==-1
break
end
        if strlength(tline)>=4
        strcheck=tline(1:4);
            if strcmp(strcheck,'  2D')

            kp=kp+1;

            %saving which parameters are being 
            %compared in this particular 2D marginal
            str_indices=tline(31:end);
            C=strsplit(str_indices);
            ind2remove=find(strlength(C)==0);
            C(ind2remove)=[];
            i=str2double(C(1));
            j=str2double(C(2));
    
            tline = fgetl(fid);
            il=il+1;
    
            tline = fgetl(fid);
            il=il+1;
    
            %saving ranges of parameter value/probability density value
            C=strsplit(tline);
            ind2remove=find(strlength(C)==0);
            C(ind2remove)=[];
            ranges(1,1,kp)=str2double(C(1));
            ranges(2,1,kp)=str2double(C(2));
            ranges(1,2,kp)=str2double(C(3));
            ranges(2,2,kp)=str2double(C(4));
    
            marg2d(1,kp) = i;
            marg2d(2,kp) = j;
            
            tline = fgetl(fid);
            il=il+1;
            %reading in 2D marginals
                for n=1:ndis
                tline = fgetl(fid);
                il=il+1;
                C=strsplit(tline);
                ind2remove=find(strlength(C)==0);
                C(ind2remove)=[];
                data(n,:,kp)=str2double(C);

                end
            end
        end
end

%final 2D marginal matrices:
data2d=data(:,:,kp-np+1:kp);
marg2d_key=marg2d(:,kp-np+1:kp); %each column of this array indicates 
                                 % which two free parameters are being 
                                 % plotted in the subplots below 
ranges2d=ranges(:,:,kp-np+1:kp); %first column of ranges(:,:,n) is the dimensional values,
                                 % second column is the scaled (non-dimensional values)

%finding minimum,maximum value for colorbar for all 2D marginal contour plots:
fmax=max(data2d,[],'all');
fmin=min(data2d,[],'all');

%plot 2-D Marginals
m=ceil(sqrt(np));  % # of rows in subplot
n=m;               % # of columns in subplot

icesheetDOF={'EIS GMSL','WAIS GMSL','EIS onset','WAIS onset','EIS duration','WAIS duration'};
figure
for p=1:np
subplot(m,n,p)
h=contourf(data2d(:,:,p));


stitle=[icesheetDOF{marg2d_key(1,p)} ' vs. ' icesheetDOF{marg2d_key(2,p)}];
xlabelstring=[icesheetDOF{marg2d_key(1,p)}];
ylabelstring=[icesheetDOF{marg2d_key(2,p)}];

xlabel(xlabelstring)
ylabel(ylabelstring)
xrange=ranges2d(:,1,p);
if marg2d_key(1,p)==1

xrange(2)=10;
end
xspacing=linspace(xrange(1),xrange(2),ndis);
xticks(linspace(0,ndis,5));
xticklabels({'0',num2str(xspacing(5)),num2str(xspacing(10)),num2str(xspacing(15)),num2str(xspacing(20))})
yrange=ranges2d(:,2,p);
yspacing=linspace(yrange(1),yrange(2),ndis);
xticks(linspace(0,ndis,5));
yticklabels({'0',num2str(yspacing(5)),num2str(yspacing(10)),num2str(yspacing(15)),num2str(yspacing(20))})
axis equal

xlim([1 ndis])
ylim([1 ndis])
colormap(flip(hot))
colorbar

end

sgtitle('Joint Posterior PDFs')

outfile='2dmarg.mat';
save(outfile,'ranges2d','ndis','marg2d_key','data2d','icesheetDOF','np')