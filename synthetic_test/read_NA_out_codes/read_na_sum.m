clear
close all

%------------------------------------------------------------------------%                                                         
%%% This script reads the summary output file from the NA Sampler      %%%    
%%% and plots the optimum model to the synthetic data and the          %%%
%%% misfit convergence over increasing NA iterations                   %%%                                      
%------------------------------------------------------------------------%

fid = fopen('na.sum');

numIce=2; % number of ice sheets considered in the MWP-1A inversion

%reading past summary text

for k=1:9
tline = fgetl(fid); %reading past summary options   
end

%get Ns, # of iterations, etc

tline = fgetl(fid); 
C = strsplit(tline);
Ns_init=str2double(C{end});

tline = fgetl(fid); 
C = strsplit(tline);
Ns=str2double(C{end});

tline = fgetl(fid); 
C = strsplit(tline);
N_iter=str2double(C{end}); % there is also a zeroth iteration

tline = fgetl(fid); 
C = strsplit(tline);
Nr=str2double(C{end});

tline = fgetl(fid); 
C = strsplit(tline);
Ntot=str2double(C{end});

%read past a few lines to next piece of info:
for k=1:8
tline = fgetl(fid);
end

%store nd, ranges, and scales
tline = fgetl(fid);
C = strsplit(tline);
nd=str2double(C{end});

for k=1:3
tline = fgetl(fid);
end

scales=zeros(nd,1);
ranges=zeros(nd,2);

for n=1:nd
tline = fgetl(fid);
C = strsplit(tline);
ranges(n,1)=str2double(C{3});
ranges(n,2)=str2double(C{4});
scales(n)=str2double(C{5});
end

for k=1:11
tline = fgetl(fid);
end

C = strsplit(tline);
min_misfit_all=str2double(C{end-1});

for k=1:2
tline = fgetl(fid);
end

C = strsplit(tline);
min_misfit_model_ix=str2double(C{end}); 

for k=1:15
tline = fgetl(fid);
end

%cycle through summary of each iteration:

for i=1:N_iter+1
tline = fgetl(fid);
C = strsplit(tline);
it_ix(i)=str2double(C{3});
tline = fgetl(fid);
C = strsplit(tline);
min_misfit(i)=str2double(C{3});

    for k=1:nd
    %read past optimum model parameters of this iteration
    tline = fgetl(fid);
    end
    
    tline = fgetl(fid); %read past empty line

end


%find optimum model:

it_to_find=floor(min_misfit_model_ix/Ns);

%read final model (optimum model of all iterations)
tline = fgetl(fid); 

model_opt=zeros(nd,1);
for n=1:nd
tline = fgetl(fid); 
C = strsplit(tline); 
model_opt(n)=str2double(C{2});
end

fname='misfit_iter.mat';

save(fname,"it_ix","min_misfit");



%% predicted data from optimum model vs actual data:
% 
%loading in synthetic data used for test inversion
observed_data=load('synthetic_data');

nMWPtsteps=37;

%% misfit values associated with each model
allout=readmatrix('predicted_data');
mfitmin=allout(end,38);
allout(:,nMWPtsteps+2:end)=[];

mfitm=allout(:,nMWPtsteps+1);

minmfit=min(mfitm);
idx=find(mfitm==minmfit);

opt_predicted(:,:)=allout(end-5:end,:);

MWPtsteps=-14.65:0.025:-13.75;

figure

for n=1:6
subplot(2,3,n)


plot(MWPtsteps,observed_data(n,1:nMWPtsteps),'bo')

hold on
plot(MWPtsteps,opt_predicted(n,1:nMWPtsteps),'rx-')
hold on

if n==1 || n==4
ylabel('sea-level change (m)')
end

set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)

xlabel('time (ka)')

if n==6
ylim([-25 25])
end


labels={'Tahiti','Barbados','Sunda Shelf','HYD (GBR)','NOG (GBR)','NW Scotland'};


title(labels{n})


end
legend({'actual data', 'predicted output of optimum model'})

sgtitle(['Site-Specific \Delta SL across MWP-1A, misfit = ' num2str(mfitmin)])
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)


%% min misfit as a function of # of iterations
figure
plot(min_misfit,'m*-')
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
ylabel('min misfit of each iteration')
xlabel('Number of NA iterations')
title('Convergence of NA sampler')

