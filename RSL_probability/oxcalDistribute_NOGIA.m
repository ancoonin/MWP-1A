function outputMatrix=oxcalDistribute_NOGIA(sampleID,noObsToGenerate,path2OxCalDist)
%------------------------------------------------------------------------%
%%% Modified from the Jan 11 2018 version from:                        %%%
%%% Hibbert, F. D., Williams, F. H., Fallon, S. J., & Rohling, E. J.   %%%
%%% Figshare https://doi.org/10.6084/m9.figshare.5890579 (2018)        %%%
%%%                                                                    %%%
%%% This function incorporates the full C-14 age distribution in the   %%%
%%% sampling of the age and RSL distributions to compute the RSL       %%%
%%% probability distribution                                           %%%
%------------------------------------------------------------------------%

%read our database data in
strIn=strcat(num2str(sampleID),'.csv');
path2file=[path2OxCalDist strIn];
customDate=readmatrix(path2file);

%here the age distribution is not necessarily Gaussian, so
%we have to consider the full C-14 age distribution 

% steps done in this code:
% 1.  use the probability distribution to form the cumulative distribution
%    cumulative distribution function is just the integral of the probability
%    distribution
% 2. check customDate for rows that are very similar
%    else get rows and rows of data that contain zeros 
%    (we do this step is because we need unique values 
%    when we interpolate the cumulative distribution
%    function to get sample times
% 3. sample the cumulative distribution function for ages in ka 

% Riemann sum: calculating the area of each rectange 
for x=2:size(customDate,1)
customDate(x,3)=(customDate(x-1,1)-customDate(x,1))*customDate(x,2);
end


% take cumulative sum of those areas to get the approx. integral of 
% the PDF from -inf to that age value

customDate(:,4)=cumsum(customDate(:,3));

for y = 2:size(customDate,1)
    if customDate(y,4)-customDate(y-1,4)<=0.000000005
        customDate(y,5)=1;
    else
        customDate(y,5)=0;
    end 
end
customDate=customDate(customDate(:,5)==0,1:4);

%note col 1 values are the times, col 4 is the CDF

sampleNos=rand(noObsToGenerate,1); % i.e. get normally distributed random values between 0 and 1

%find the values within the age distribution 
% that these points correspond to:

%i.e. sampling the cumulative distribution function for times in ka 

outputMatrix=(((interp1(customDate(:,4),customDate(:,1),...
    sampleNos)))/1000); 

end



