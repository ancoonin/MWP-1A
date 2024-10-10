function [elevation_upliftCorrected,elevationError_upliftCorrected ] =...
    upliftCorrectElevation_NOGIA...
    ( coralData,ageApplication,rowIdentifierCorals,rowIdentifierAge)

%------------------------------------------------------------------------%
%%% Modified from the Jan 11 2018 version from:                        %%%
%%% Hibbert, F. D., Williams, F. H., Fallon, S. J., & Rohling, E. J.   %%%
%%% Figshare https://doi.org/10.6084/m9.figshare.5890579 (2018)        %%%
%%%                                                                    %%%
%%% This function corrects elevations for corals given the             %%%
%%% corresponding age sampled from the age distribution and            %%%
%%% propagates the associated error into the RSL uncertainty           %%%
%------------------------------------------------------------------------%

    if coralData(rowIdentifierCorals,3)~=0
        %i.e. where uplift has a non-zero value
        elevation_upliftCorrected=coralData(rowIdentifierCorals,5)-...
            (coralData(rowIdentifierCorals,3) *...
            ageApplication(rowIdentifierAge,1));
        errorTerm1=coralData(rowIdentifierCorals,6)^2;
        errorTerm2=(coralData(rowIdentifierCorals,3)*...
            ageApplication(rowIdentifierAge,1))*sqrt...
            (((coralData(rowIdentifierCorals,4)/coralData...
            (rowIdentifierCorals,3))^2)+...
            (((coralData(rowIdentifierCorals,12)/2)/...
            ageApplication(rowIdentifierAge,1))^2));
        elevationError_upliftCorrected=sqrt((errorTerm1) + (errorTerm2^2));
    else
        %can keep same equation for zedtcp, but not for error
elevation_upliftCorrected=coralData(rowIdentifierCorals,5)-...
    (coralData(rowIdentifierCorals,3) *...
    ageApplication(rowIdentifierAge,1));
elevationError_upliftCorrected=coralData(rowIdentifierCorals,6);
    end 
end


