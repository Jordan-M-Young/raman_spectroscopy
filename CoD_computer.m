function [CoD_R2,ssss] = CoD_computer(extracted_data,raw_data)
% This function requires an extracted_data array and a raw_data array
% make sure that these arrays are the same size!

% CoD_R2 and CoD_R2a will hold R^2 values that are computed from different
% formulae, but should be equivalent 
CoD_R2 = zeros(1,45000);
CoD_R2 = transpose(CoD_R2);
CoD_R2a = zeros(1,45000);
CoD_R2a = transpose(CoD_R2a);
[rows columns] = size(raw_data);


elements = rows*columns;
raw_shift = zeros(1,columns);
extracted_shift = zeros(1,columns);

for i = 1:rows
    raw_shift = raw_data((1*i):rows:((elements-rows) +(1*i)));
    extracted_shift = extracted_data((1*i):rows:((elements-rows) +(1*i)));
    
    %CoD or R^2 calculation based on:
    % https://newonlinecourses.science.psu.edu/stat501/node/255/
  
    y_hat = extracted_shift;
    y = raw_shift;
    y_bar = mean(y_hat);
    
    SSTO = (y - y_bar).^2;
    SSTO = sum(SSTO);
    
    SSR = (y_hat - y_bar).^2;
    SSR = sum(SSR);
 
   
    SSE = (y - y_hat).^2;
    SSE = sum(SSE);
    
    CoD_R2(i) = 1 - (SSE / SSTO);
    CoD_R2a(i) = (SSR / SSTO);
    
end
suc = zeros(1,45000);

%checks to see if the independent calculations 
% of a samples R^2 value are equivalent
for i = 1:size(CoD_R2)
    if CoD_R2(i) - CoD_R2a(i) < 0.1;
        suc(i) = 1;
    end
end

    
CoD_R2
suc = sum(suc)
