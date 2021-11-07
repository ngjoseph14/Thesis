
close all;
clear all;

% ------------------------------------------------------------------------

% These two lines of code ensure the graphs have a white bakcground 
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% Inputs                                
Window_length = 256;        % Sets the window length of the segments of data

% ------------------------------------------------------------------------

% Load Dataset
Dataset = load('Altera_CWLFM_Dataset.mat');
HRR_profiles = Dataset.HRR_profiles;                    % Obtains the 2D matrix of the data (with the rnage profiles and range bins)
ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;    % Obtains the PRF for the dataset
Range_axis = Dataset.Range_axis;                        % Obtains the array for the range axis

% ------------------------------------------------------------------------

% Plot All HRR profiles
NumProfiles_axis = 1:size(HRR_profiles,1);      % Obtains the range profile array
figure(); 
imagesc(Range_axis, NumProfiles_axis, 20*log10(abs(HRR_profiles))); 
colorbar;
xlabel('Range (m)');        % X-axis label in range
ylabel('Profile number');   % Y-axis label in range profiles
title('Unaligned range profiles');      % Title for the plot
axis xy;
colormap('jet');

% ------------------------------------------------------------------------

overlap = 0.5*Window_length;                    % Computes overlap amount/factor (data points reused in consecutive windows)
hop_length = Window_length - overlap;           % Computes hop length (data points between window segments)
Image_count = floor((size(HRR_profiles,1) - overlap)/hop_length);       % Computes the number of images (sets) that the iteration will run through and generate
start_profile = 1;  

for o = 1:Image_count
    Subset_HRR_profiles = HRR_profiles((start_profile):(Window_length + start_profile - 1),:);      % Subdivides the data into segments to be run through the ISAR processor

    [ref_profile, n_order] = find_Reference_profile_and_n_order(Subset_HRR_profiles, Window_length, o);     % Finds the best reference profile and polynomial order number to produce the best image contrast value
    
    Aligned_HRR_profiles = Haywood_Range_Alignment(Subset_HRR_profiles, ref_profile, n_order);              % Applies the Haywood range alignment to the segmented data
    ref_bin = DSA_Phase_Adjustment_selected_scatterer_plot(Aligned_HRR_profiles);                           % Phase compensation is applied to the range aligned data to find the reference range bins that meet the first two criteria 
    
    [top_ref_bin, Top_MSE] = find_lowest_MSE(Aligned_HRR_profiles, ref_bin);                % The five range bins that results in teh lowest MSE value are obtained
    best_ref_bin = find_best_ref_bin(Aligned_HRR_profiles, top_ref_bin, Window_length);     % The range bin that results in the best IC value is chosen as the reference range bin  
    fprintf('Reference range bin for Image %i with DSA is: %f \n', o, best_ref_bin);        % Displays the chosen reference range bin
    plot_phase(Aligned_HRR_profiles, best_ref_bin);                                         % Plots the phase of the scatterer in the chosen reference range bin
    
    PhaseComp_HRR_profiles = DSA_Phase_Adjustment(Aligned_HRR_profiles, best_ref_bin);      % Phase compensation is applied to the range aligned data with the reference range bin
    
    Windowing = repmat(hamming(Window_length),1,size(Aligned_HRR_profiles,2));                      % A Hamming windowing function (same size as the segmented data) is created
    Windowed_HRR_profiles = PhaseComp_HRR_profiles.*Windowing;                                      % The Hamming windowing function is applied to the phase compensated data 
    Windowed_noDSA_HRR_profiles = Aligned_HRR_profiles.*Windowing;                                  % The Hamming windowing function is applied to the range aligned data
    
    DopplerAnalysis_HRR_profiles = fftshift(fft(Windowed_HRR_profiles,[],1), 1);                    % Doppler processing is applied to the windowed phase compensated data (to obtain final ISAR images)
    DopplerAnalysis_noDSA_HRR_profiles = fftshift(fft(Windowed_noDSA_HRR_profiles,[],1), 1);        % Doppler processing is applied to the windowed range aligned data (to obtain the ISAR images before autofocus)
    
    calculate_contrast = calculate_image_contrast(DopplerAnalysis_HRR_profiles);                    % Calculates the image contrast for the ISAR images after autofocus
    fprintf('Image Contrast for Image %i with DSA is: %f \n', o, calculate_contrast);               % Displays the calculated image contrast value from above
    
    calculate_contrast_noDSA = calculate_image_contrast(DopplerAnalysis_noDSA_HRR_profiles);        % Calculates the image contrast for the ISAR images before autofocus
    fprintf('Image Contrast for Image %i without DSA is: %f \n', o, calculate_contrast_noDSA);      % Displays the calculated image contrast value from above
    
    % Plot unaligned HRR profiles for the segment of data
    unalignedNumProfiles_axis = 1:size(Subset_HRR_profiles,1);      % Obtains the range profile array (y-axis)
    figure(); 
    imagesc(Range_axis, unalignedNumProfiles_axis, 20*log10(abs(Subset_HRR_profiles))); 
    colorbar;
    xlabel('Range (m)');                % X-axis label in range
    ylabel('Profile number');           % Y-axis label in range profiles
    title('Unaligned range profiles');  % Title for the plot
    axis xy;
    colormap('jet');
    
    % Plot aligned HRR profiles for the segment of data
    alignedNumProfiles_axis = 1:size(Aligned_HRR_profiles,1);       % Obtains the range profile array (y-axis)
    figure(); 
    imagesc(Range_axis, alignedNumProfiles_axis, 20*log10(abs(Aligned_HRR_profiles))); 
    colorbar;
    xlabel('Range (m)');                % X-axis label in range
    ylabel('Profile number');           % Y-axis label in range profiles
    title('Aligned range profiles');    % Title for the plot
    axis xy;
    colormap('jet');
    
    % Computation to convert the y-axis from range profiles to Doppler frequency
    n = size(HRR_profiles,1);
    y_axis = ProfileRepetitionFreq*(Window_length/2:-1:(-(Window_length/2-1)))/Window_length;       % Y-axis in terms of Doppler frequency (for the ISAR images)
    
    % This code is used to make the ISAR image results more clear to see
    peakVal_noDSA = max(max(abs(DopplerAnalysis_noDSA_HRR_profiles)));
    Normalise_ISAR_image_noDSA_abs = abs(DopplerAnalysis_noDSA_HRR_profiles)/peakVal_noDSA;
    clims_noDSA = [-50 0];                                                  % Only plot values from -50 dB to 0 dB
    
    % Plot the ISAR images of the segment of data before autofocus
    figure();
    imagesc(Range_axis, y_axis, 20*log10(abs(Normalise_ISAR_image_noDSA_abs)), clims_noDSA); 
    colorbar;
    xlabel('Range (m)');                % X-axis label in range
    ylabel('Doppler frequency (Hz)');   % Y-axis label in Doppler frequency
    title(['ISAR Image ', num2str(o) ,' without DSA']);     % Title for the plot
    axis xy;
    colormap('jet');
    
    % This code is used to make the ISAR image results more clear to see
    peakVal = max(max(abs(DopplerAnalysis_HRR_profiles)));
    Normalise_ISAR_image_abs = abs(DopplerAnalysis_HRR_profiles)/peakVal;
    clims = [-50 0];                                                        % Only plot values from -50 dB to 0 dB
    
    % Plot the ISAR images of the segment of data before autofocus
    figure();
    imagesc(Range_axis, y_axis, 20*log10(abs(Normalise_ISAR_image_abs)), clims); 
    colorbar;
    xlabel('Range (m)');                % X-axis label in range
    ylabel('Doppler frequency (Hz)');   % Y-axis label in Doppler frequency
    title(['ISAR Image ', num2str(o) ,' with DSA']);    % Title for the plot
    axis xy;
    colormap('jet');
    
    start_profile = start_profile + hop_length;         % Updates the starting profile for the next window segment
end

% ------------------------------------------------------------------------

function [Ref_HRR_profile, n_order] = find_Reference_profile_and_n_order(HRR_profiles, Window_length, image_count)
% This function finds the optimal reference range profile and polynomial order to produce the best IC value

N = size(HRR_profiles,1);       % Obtains the number of the range profiles
contrast_max = 0;

for f = 1:N
    for h = 1:1:3
        Aligned_HRR_profiles = Haywood_Range_Alignment(HRR_profiles, f, h);                 % Aligns range profiles for each range profile value and polynomial order value     
    
        % The process of phase compensation is then applied
        ref_bin = DSA_Phase_Adjustment_selected_scatterer(Aligned_HRR_profiles);            % Finds the range bins that meet the first two criteria       
        [top_ref_bin, Top_MSE] = find_lowest_MSE(Aligned_HRR_profiles, ref_bin);            % The five range bins that results in teh lowest MSE value are obtained
        best_ref_bin = find_best_ref_bin(Aligned_HRR_profiles, top_ref_bin, Window_length); % The range bin that results in the best IC value is chosen as the reference range bin
        PhaseComp_HRR_profiles = DSA_Phase_Adjustment(Aligned_HRR_profiles, best_ref_bin);  % Phase compensation is applied to the range aligned data with the reference range bin
    
        Windowing = repmat(hamming(Window_length),1,size(Aligned_HRR_profiles,2));          % A Hamming windowing function (same size as the segmented data) is created
        Windowed_HRR_profiles = PhaseComp_HRR_profiles.*Windowing;                          % The Hamming windowing function is applied to the phase compensated data
    
        DopplerAnalysis_HRR_profiles = fftshift(fft(Windowed_HRR_profiles,[],1), 1);        % Doppler processing is applied to the windowed phase compensated data
    
        calculate_contrast = calculate_image_contrast(DopplerAnalysis_HRR_profiles);        % Computes the IC value for the iteration
        if calculate_contrast > contrast_max                                                % Compares to the current max IC value
            contrast_max = calculate_contrast;                                              % If the new IC value is greater, it becomes the max IC value
            Ref_HRR_profile = f;                                                            % Updates the reference range profile variable
            n_order = h;                                                                    % Updates the reference range profile variable
        end
    end
end

fprintf('Reference Profile for Image %i in ISAR Processor is: %f \n', image_count, Ref_HRR_profile);    % Displays the selected reference range profile
fprintf('Order of polynomial for Image %i in ISAR Processor is: %f \n', image_count, n_order);          % Displays the selected polynomial order

end

function Aligned_HRR_profiles = Haywood_Range_Alignment(HRR_profiles, ref_profile, n_order)
% This function is implemented to perform the Haywood range alignment technique

N = size(HRR_profiles,1);   % Obtains the number of the range profiles     
M = size(HRR_profiles,2);   % Obtains the number of the range bins         
i = (0:1:(M-1));                 
Aligned_HRR_profiles = zeros(N,M);      % Creates an array of size NxM for the range aligned data

for a = 1:N
    [r,lags] = xcorr(abs(HRR_profiles(ref_profile,:)),abs(HRR_profiles(a,:)));            % Returns cross correlation and lags values for each computation between reference profile (first profile) and the chosen profile 
    peak = max(r);                                                              % Finds max value from the cross correlation values
    index_peak = find(r == peak);                                               % Finds index of the max value
    
    delay = lags(index_peak);                                                   % Obtains the delay value of the peak
    delayVector(a) = delay;                                                     % Saves the delay values per profile in an array (this is before the curve fitting of the polynomial)
    lin_delay(a) = a;                                                           % Creates the linspace for the delay vector
end

poly_coeff = polyfit(lin_delay,delayVector,n_order);                            % Finds the coefficients for the n order polynomial
smoothed_delays = polyval(poly_coeff, lin_delay).';                             % Calculates the new smoothed delay vector
   
for b = 1:N
    shift_amount = exp(-1i*(2*pi*smoothed_delays(b).*i/M));                     % Required shift needed for the current profile selected
    Aligned_HRR_profiles(b,:) = ifft(shift_amount.*fft(HRR_profiles(b,:)));     % Apply the shift to the current profile selected in the frequency domain
end

end

function ref_bin = DSA_Phase_Adjustment_selected_scatterer(Aligned_HRR_profiles)
% This function finds the rnage bins that meet the first two criterion

N = size(Aligned_HRR_profiles,1);           % Obtains the number of the range profiles
T = size(Aligned_HRR_profiles,2);           % Obtains the number of the range bins 

Var_array = var(abs(Aligned_HRR_profiles));                                 % Creates an array with all the variances of each range bin
Sort_var = sort(Var_array);                                                 % Sorts the variances from min value to max value

Power_vector = sum(abs(Aligned_HRR_profiles).^2);                           % Creates an array with the sum of the total power of each range bin over all HRRPs
Total_power = sum(Power_vector)/N;                                          % The value for the total power over all range bins and all HRRPs                                         

min_var = 0;
j = 1;

for k = 1:size(Aligned_HRR_profiles,2)
    min_var = Sort_var(k);                                                  % Selects a variance value starting from the minimum value to maximum value
    index_bin = find(Var_array == min_var);                                 % Finds the index of the range bin for the selected variance
    Power_bin = Power_vector(index_bin);                                    % Finds the total power for that range bin over all HRRPs
    if Power_bin > Total_power                                              % Checks whether the selected range bin meets the dominant scatterer criterion
        ref_bin(j) = index_bin;                                             % Adds the range bin to an array
        j = j + 1;                                                          % Continues to next iteration
    end
end

end

function ref_bin = DSA_Phase_Adjustment_selected_scatterer_plot(Aligned_HRR_profiles)

N = size(Aligned_HRR_profiles,1);       % Obtains the number of the range profiles
T = size(Aligned_HRR_profiles,2);       % Obtains the number of the range bins 

Var_array = var(abs(Aligned_HRR_profiles));                                 % Creates an array with all the variances of each range bin
Sort_var = sort(Var_array);                                                 % Sorts the variances from min value to max value

Power_vector = sum(abs(Aligned_HRR_profiles).^2);                           % Creates an array with the sum of the total power of each range bin over all HRRPs
Total_power = sum(Power_vector)/N;                                          % The value for the total power over all range bins and all HRRPs                                         

min_var = 0;
j = 1;

for k = 1:size(Aligned_HRR_profiles,2)
    min_var = Sort_var(k);                                                  % Selects a variance value starting from the minimum value to maximum value
    index_bin = find(Var_array == min_var);                                 % Finds the index of the range bin for the selected variance
    Power_bin = Power_vector(index_bin);                                    % Finds the total power for that range bin over all HRRPs
    if Power_bin > Total_power                                              % Checks whether the selected range bin meets the dominant scatterer criterion
        ref_bin(j) = index_bin;                                             % Selects the range bin as the reference bin if criterion is met
        j = j + 1;                                                          % Continues to next iteration of finding the reference bin if criterion is not met 
    end
end

lin_space = 1:1:T;
Power_threshold = repmat(Total_power,T);

% Plots the amplitude variance for all of the range bins
figure();
plot(lin_space, Var_array);
hold on
plot(ref_bin, Var_array(ref_bin),'r.','MarkerSize',15);
xlabel('Range Bins');       % X-axis label in range bins
ylabel('Variance');         % Y-axis label in amplitude variance
legend('Range bin Variance');

% Plots the amplitude power for all of the range bins
figure();
plot(lin_space, Power_vector);
hold on
plot(lin_space, Power_threshold,'m--');
hold on
plot(ref_bin, Power_vector(ref_bin),'r.','MarkerSize',15);
xlabel('Range Bins');       % X-axis label in range bins
ylabel('Amplitude Power');  % Y-axis label in amplitude power
legend('Range bin Power','Threshold Power');

end

function PhaseComp_HRR_profiles = DSA_Phase_Adjustment(Aligned_HRR_profiles, ref_bin_input)
% This function computes the phase compensation

N = size(Aligned_HRR_profiles,1);       % Obtains the number of the range profiles
Dominant_scatterer = Aligned_HRR_profiles(:,ref_bin_input);     % Obtains the dominant scatterer vector needed to compute the compensation vector
Complex_Conjugate = conj(Dominant_scatterer);       % Computes the compensation vector by taking the complex conjugate of the dominant scatterer

for m = 1:N                                    
    PhaseComp_HRR_profiles(m,:) = Aligned_HRR_profiles(m,:) * Complex_Conjugate(m);   % Apply the compensation vector to the range profiles in the segmented data
end

end

function calculate_contrast = calculate_image_contrast(ISAR_image_complex_linear)
% This function computes the image contrast for the ISAR images 

I_temp = abs(ISAR_image_complex_linear).^2;
mean_mean_I_temp = mean( mean(I_temp ));


NumeratorPart1 = I_temp - mean_mean_I_temp;
NumberatorFull = sqrt(mean(mean(NumeratorPart1.^2)));

Image_contrast = NumberatorFull./mean_mean_I_temp;
calculate_contrast = Image_contrast;

end

function [top_ref_bin, MSE_array] = find_lowest_MSE(Aligned_HRR_profiles, ref_bin_input)

    N = size(Aligned_HRR_profiles,1);       % Obtains the number of the range profiles

    for h = 1:length(ref_bin_input)
        phase_vector = angle(Aligned_HRR_profiles(:,ref_bin_input(h)));     % Obtains the phase of the dominant scatterer for each iteration
        linspace_phase = 1:N;
        phase_vector_unwrapped = unwrap(phase_vector);                      % Unwraps the phase of the domianant scatterer for each iteration (to allow for a linear regression line to be compared to it)
    
        poly_coeff = polyfit(linspace_phase, phase_vector_unwrapped, 1);    % Finds the coefficients for the first order polynomial (linear regression line) for each iteration                                   
        linear_progression = polyval(poly_coeff, linspace_phase).';         % Obtains the linear regression line (vector) over the linspace for each iteration  
    
        MSE_array(h) = immse(linear_progression ,phase_vector_unwrapped);   % Computes the MSE value for each iteration
    end
   
    sort_MSE = sort(MSE_array);     % Sorts the MSE array from min to max

    top_MSE = sort_MSE(1:5);        % Obtains the 5 lowest MSE values
    for q = 1:length(top_MSE)
        min_MSE = sort_MSE(q);      % Selects the minimum MSE value
        index_MSE = find(MSE_array == min_MSE);     % Obtains the index of this value from the original MSE array
        top_ref_bin(q) = ref_bin_input(index_MSE);  % Stores the range bin that corresponds to that MSE value (index)
    end
end

function plot_phase(Aligned_HRR_profiles, top_ref_bin)
% This function plots the phase of the scatterer in the selected reference range bin

    N = size(Aligned_HRR_profiles,1);       % Obtains the number of the range profiles
    
    phase_vector = angle(Aligned_HRR_profiles(:,top_ref_bin));  % Obtains the phase of the dominant scatterer 
    linspace_phase = 1:N;
    phase_vector_unwrapped = unwrap(phase_vector);              % Unwraps the phase of the domianant scatterer (to allow for a linear regression line to be compared to it)
    
    poly_coeff = polyfit(linspace_phase, phase_vector_unwrapped, 1);    % Finds the coefficients for the first order polynomial (linear regression line)                                  
    linear_progression = polyval(poly_coeff, linspace_phase).';         % Obtains the linear regression line (vector) over the linspace 
    
    % Plots the phase of the dominant scatterer versus the linear regression line
    figure();
    plot(linspace_phase, phase_vector_unwrapped,'r--','LineWidth',2);
    hold on;
    plot(linspace_phase, linear_progression,'b:','LineWidth',2);
    xlabel('Range Profile Number');     % X-axis label in range profiles
    ylabel('Phase');                    % Y-axis label in phase
    legend('Phase of scatterer','linear curve fit');
    
end

function best_ref_bin = find_best_ref_bin(Aligned_HRR_profiles, top_ref_bin, Window_length)
% This function finds the rnage bin that produces the best IC value

 max_constrast = 0;
    
    for r = 1:length(top_ref_bin)
        PhaseComp_HRR_profiles = DSA_Phase_Adjustment(Aligned_HRR_profiles, top_ref_bin(r));    % Phase compensation is applied to the range aligned data with the reference range bin
    
        Windowing = repmat(hamming(Window_length),1,size(Aligned_HRR_profiles,2));          % A Hamming windowing function (same size as the segmented data) is created
        Windowed_HRR_profiles = PhaseComp_HRR_profiles.*Windowing;                          % The Hamming windowing function is applied to the phase compensated data
    
        DopplerAnalysis_HRR_profiles = fftshift(fft(Windowed_HRR_profiles,[],1), 1);        % Doppler processing is applied to the windowed phase compensated data
        calculate_contrast = calculate_image_contrast(DopplerAnalysis_HRR_profiles);        % Computes the IC value
        
        if calculate_contrast > max_constrast                                               % Checks to see if the IC value is greater than the current max IC value
            best_ref_bin = top_ref_bin(r);                                                  % Updates the reference range bin if met
            max_constrast = calculate_contrast;                                             % Updates the max contrast value if met as well
        end
    end
end