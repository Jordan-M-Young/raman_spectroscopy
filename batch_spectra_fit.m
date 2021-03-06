%loads your raw hyperspectral data-cube .mat file; enter directory/filename here
load('SmoothedMaps/Raw/ALH77012(1)(SG).mat') 

%loads your fit hyperspectral data-cube .mat file; enter directory/filename here
load('SmoothedMaps/Extracted/ALH77012(1)(SG)EXTRACTED.mat') 

%accesses spectral array from raw data structure
raw_data = ALH770121SG.data;

%accesses spectral array from fit data structure
data =  ALH770121SGEXTRACTED.data;

[row,col] = size(data)

% Wavenumber values for spectral channels of interest
wavenumber = [9.04E+02	9.08E+02	9.12E+02	9.17E+02	9.21E+02	9.25E+02	9.30E+02	9.34E+02	9.38E+02	9.43E+02	9.47E+02	9.52E+02	9.56E+02	9.60E+02	9.65E+02	9.69E+02	9.73E+02	9.78E+02	9.82E+02	9.86E+02	9.90E+02	9.95E+02	9.99E+02	1.00E+03	1.01E+03	1.01E+03	1.02E+03	1.02E+03	1.02E+03	1.03E+03	1.03E+03	1.04E+03	1.04E+03	1.05E+03	1.05E+03	1.06E+03	1.06E+03	1.06E+03	1.07E+03	1.07E+03	1.08E+03	1.08E+03	1.08E+03	1.09E+03	1.09E+03	1.10E+03	1.10E+03	1.11E+03	1.11E+03	1.11E+03	1.12E+03	1.12E+03	1.13E+03	1.13E+03	1.14E+03	1.14E+03	1.14E+03	1.15E+03	1.15E+03	1.16E+03	1.16E+03	1.17E+03	1.17E+03	1.17E+03	1.18E+03	1.18E+03	1.19E+03	1.19E+03	1.20E+03	1.20E+03	1.20E+03	1.21E+03	1.21E+03	1.22E+03	1.22E+03	1.22E+03	1.23E+03	1.23E+03	1.24E+03	1.24E+03	1.25E+03	1.25E+03	1.25E+03	1.26E+03	1.26E+03	1.27E+03	1.27E+03	1.27E+03	1.28E+03	1.28E+03	1.29E+03	1.29E+03	1.30E+03	1.30E+03	1.30E+03	1.31E+03	1.31E+03	1.32E+03	1.32E+03	1.32E+03	1.33E+03	1.33E+03	1.34E+03	1.34E+03	1.35E+03	1.35E+03	1.35E+03	1.36E+03	1.36E+03	1.37E+03	1.37E+03	1.37E+03	1.38E+03	1.38E+03	1.39E+03	1.39E+03	1.40E+03	1.40E+03	1.40E+03	1.41E+03	1.41E+03	1.42E+03	1.42E+03	1.42E+03	1.43E+03	1.43E+03	1.44E+03	1.44E+03	1.44E+03	1.45E+03	1.45E+03	1.46E+03	1.46E+03	1.46E+03	1.47E+03	1.47E+03	1.48E+03	1.48E+03	1.49E+03	1.49E+03	1.49E+03	1.50E+03	1.50E+03	1.51E+03	1.51E+03	1.51E+03	1.52E+03	1.52E+03	1.53E+03	1.53E+03	1.53E+03	1.54E+03	1.54E+03	1.55E+03	1.55E+03	1.55E+03	1.56E+03	1.56E+03	1.57E+03	1.57E+03	1.57E+03	1.58E+03	1.58E+03	1.59E+03	1.59E+03	1.59E+03	1.60E+03	1.60E+03	1.61E+03	1.61E+03	1.61E+03	1.62E+03	1.62E+03	1.63E+03	1.63E+03	1.63E+03	1.64E+03	1.64E+03	1.65E+03	1.65E+03	1.65E+03	1.66E+03	1.66E+03	1.67E+03	1.67E+03	1.67E+03	1.68E+03	1.68E+03	1.69E+03	1.69E+03	1.69E+03	1.70E+03	1.70E+03	1.71E+03	1.71E+03	1.71E+03	1.72E+03	1.72E+03	1.73E+03	1.73E+03	1.73E+03	1.74E+03	1.74E+03	1.75E+03	1.75E+03	1.75E+03	1.76E+03	1.76E+03	1.76E+03	1.77E+03	1.77E+03	1.78E+03	1.78E+03	1.78E+03	1.79E+03	1.79E+03	1.80E+03	1.80E+03	1.80E+03	1.81E+03	1.81E+03	1.82E+03	1.82E+03	1.82E+03	1.83E+03	1.83E+03	1.84E+03	1.84E+03	1.84E+03	1.85E+03	1.85E+03	1.85E+03	1.86E+03	1.86E+03	1.87E+03	1.87E+03	1.87E+03	1.88E+03	1.88E+03	1.89E+03	1.89E+03	1.89E+03	1.90E+03	1.90E+03	1.91E+03	1.91E+03	1.91E+03	1.92E+03	1.92E+03	1.92E+03	1.93E+03	1.93E+03	1.94E+03	1.94E+03	1.94E+03	1.95E+03	1.95E+03	1.96E+03	1.96E+03	1.96E+03	1.97E+03	1.97E+03	1.97E+03	1.98E+03	1.98E+03	1.99E+03	1.99E+03	1.99E+03	2.00E+03	2.00E+03	2.00E+03	2.01E+03	2.01E+03	2.02E+03];

% number of raman bands you wish to fit
numbands = 2


% array of zeros that will eventually hold your spectral parameters
% array has dimensions numbands*5 x row, each band has 5 attributes 
% and row is the number of pixel/spectra pairs in your hyperspectral image
carbon_spectra_params = zeros(numbands*5,row);

%array of zeros that will eventually hold your spectral data
carbon_spectra = zeros(col,row);

%Calculates the coefficient of determination (R^2) value
% between pairs of fit/raw spectra between data and raw_data
%CoD_R2 is an array of R^2 values and suc is an integer value
% of succesful computations of meaningful R^2 values
[CoD_R2,suc] = CoD_computer(data,raw_data);

% an array of zeros that will hold the new index numbers
% of relevant spectra
counting = zeros(1,row);

%an array of zeros that will hold the maximum values of
%each relevant spectra for normalization purposes
Normalize_vals = zeros(1,row);

%sets a counter for reindexing relevant pixels/spectra into carbon_spectra_params
% and carbon_spectra arrays
counter = 1;

%This block does several things:
for i = 1:row
    % discriminates against fit spectra based on their R^2 value. Set it to what you want!
    if CoD_R2(i) > 0.95
        
        % pulls spectrum intensity values for sample 'i' 
        intensity = data((1*i):row:((row*col) - row +(1*i)));
        intensity = double(intensity);
        
        % signal is a 2 x num spectral channels array with wavenumbers as row 1 and intensities as row 2
        % this is our spectrum
        signal = [wavenumber;intensity];
        
        %fits spectrum; FitResults are the spectral parameters you want
        % see https://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm for more information
        % Info about the author of this function:
        % This code is part of "A Pragmatic Introduction to Signal Processing", created and maintained by 
        %Prof. Tom O'Haver , Department of Chemistry and Biochemistry, The University of Maryland at College Park. 
        %Comments, suggestions and questions should be directed to Prof. O'Haver at toh@umd.edu. 
        [FitResults,FitError] = peakfit(signal,1400,1100,2,2,1);
        
        % screens spectral with peak positions (cm-1) that are outside realistic bounds that you determine
        % add more if you have more than 2 bands
        
        %    band1 lower bound       band1 upper bound       band2 lower bound       band2 upper bound
        if  FitResults(2) <= 1300 | FitResults(2) >= 1400 | FitResults(7) <= 1550 | FitResults(7) >= 1620 
        else   
            FitResults = transpose([FitResults(1:2:9),FitResults(2:2:10)]);
            
            % adds the extracted spectral parameters for sample 'i' to carbon_spectra_params
            carbon_spectra_params(1 + ((counter -1)*10):10+((counter -1)*10)) = FitResults;
            
            % adds the fit spectral intensities for sample 'i' to carbon_spectra
            carbon_spectra(1 + ((counter -1)*col):col+((counter -1)*col)) = transpose(y);
            
            % adds index to new set of indices
            counting(counter) = i ;
            
            % adds maximum intensity of sample 'i' to Normalize_vals
            Normalize_vals(counter) = max(intensity);
            
            counter = counter + 1
            
        end
    else       
    end
end

%adds old indices to new matrix
carbon_spectra_params = [carbon_spectra_params;counting];


%writes spectral parameters to a .csv file of name you choose.
csvwrite('ALH770121_95_SG_params.csv',carbon_spectra_params)

%writes spectral intensities to a .csv file of name you choose.
csvwrite('ALH770121_95_SG_spectra.csv',carbon_spectra)
