function outputArray=coralSimulations(sectionedCorals,taxaCount,...
    taxaList,timeInstances,rslInstances,buildCase,path2CoralEmpiricalDepthDist)
%------------------------------------------------------------------------%
%%% Modified from the Jan 11 2018 version from:                        %%%
%%% Hibbert, F. D., Williams, F. H., Fallon, S. J., & Rohling, E. J.   %%%
%%% Figshare https://doi.org/10.6084/m9.figshare.5890579 (2018)        %%%
%%%                                                                    %%%
%%% This function incorporates the empirical depth distribution from   %%%
%%% of modern coral analogues for coral data                           %%%
%------------------------------------------------------------------------%

    %how many of the similar taxa
    NoTaxaCorals=size(sectionedCorals,1);


    
    switch(buildCase)
        case 1
            path2file=[path2CoralEmpiricalDepthDist 'Global/'];
            %cd('Global Depth Distribution Extracts Locations/');
        case 2
            path2file=[path2CoralEmpiricalDepthDist 'Regional/'];
            %cd('/Regional Depth Distribution Extracts Locations/');
    end

    filename=[path2file num2str(taxaList(taxaCount,1)) '.csv'];

    %exactDistributionImport=-1*(csvread(filename));
     exactDistributionImport=-1*(readmatrix(filename));

    coralFit = fitdist(exactDistributionImport, 'kernel');
    outputArray=random(coralFit,NoTaxaCorals,timeInstances,rslInstances);
    %this needs normalising
    myMedian=median(reshape(outputArray,...
        NoTaxaCorals*rslInstances*timeInstances,1));
    outputArray=outputArray-myMedian;
end

