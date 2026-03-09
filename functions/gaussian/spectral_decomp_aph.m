function [amps, compspec, sumspec] = spectral_decomp_aph(wl,aph,wlunc,apunc,acs)

% AC-S "un-smoothing" and spectral decomposition method
%
% Ron Zaneveld, WET Labs, Inc., 2005
% Ali Chase, University of Maine, 2014
% 
% Vectorized: Guillaume Bourdin, University of Maine, 2020
%
% See the following publication for details on the method:
% Chase, A., et al., Decomposition of in situ particulate absorption
% spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
%
% INPUTS:
% wl     -  wavelengths associated with measured particulate absorption
%           spectra
% ap     -  measured particulate absorption spectra
% wlunc  -  wavelengths associated with uncertainty of measured
%           particulate absorption spectra
% apunc  -  uncertainty of measured particulate absorption spectra
% acs    -  1 (for AC-S data) or 0 (for a different instrument). This code
%           was designed for use with particulate absorption spectra measured
%           using a WETLabs AC-S instrument, and thus it applies a spectral
%           un-smoothing that is specific to that instrument. To use with
%           other types of absorption data, input "0" and the correction
%           applied for AC-S will be bypassed.
%
% The uncertainty values we use are the standard deviation of one-minute
% binned particulate absorption spectra. If the uncertainty is unknown,
% then all wavelenghts will be treated equally, i.e. no spectral weighting
% will be applied. In this case populate apunc with -999 values.
%
% OUTPUTS:
% amps      -  the amplitudes of 12 Gaussian function peaks and one non-algal
%              particle (NAP) function
% compspec  -  the component spectra after spectral decomposition; the
%              Gaussian and NAP functions multiplied by the 'amps'
%              values
% sumspec   -  the sum of the component spectra. will be simliar to the
%              original measured spectrum, but the fit may not be as good in
%              parts of the spectrum with higher uncertainty.
%
% The amplitudes of the Gaussian functions represent the relative amounts
% of different phytoplankton pigments and can be compared to HPLC pigment
% concentrations to evaluate the method. See Chase et al. (2014) for
% detail on such evaluation (reference above).

% if any(min(apunc,[],2) < 0)
%     apunc(min(apunc,[],2) < 0, :) = ones(max(length(ap)), size(apunc(min(apunc,[],2) < 0, :), 1));
% end
% 
% if size(apunc,1) < size(apunc,2)
%     apunc = apunc';
% end

if acs == 1

    fprintf('Absorption spectra correction ... ')

    % Set up filter factors at every 0.1 nm from 1 to 799 nm, with center
    % wavelength at centwavel (i.e. at the data wavelengths)
    wavelength = .1:.1:799; % Thus index of 1 nm = 10; 356 nm= 3560;
    clear filtfunc
    SIG1 = (-9.845*10^-8.*wl.^3  +1.639*10^-4*wl.^2- 7.849*10^-2*wl + 25.24)/2.3547 ;
    for i = 1:max(size(wl))
        for jkl = 1:max(size(wavelength))
            filtfunc(jkl,i) = (1/(sqrt(2*pi)*SIG1(i)))*exp(-0.5*((wavelength(jkl)-wl(i))/SIG1(i)).^2); % First term normalizes area under the curve to 1.
        end
    end

    % Convolve the measurement with the fiter factors add the difference to
    % the measured spectrum to get the first corrected spectrum.
    % This is the corrected absorption spectrum "abscorr".

    minwavel = min(wl);
    maxwavel = max(wl);

    % Define xixi
    xixi = minwavel:.1:maxwavel;

    % Initialize yiyi with NaNs
    yiyi = NaN(size(ap,1), length(xixi));

    % Perform interpolation row-wise
    for i = 1:size(ap,1) % Loop over each station
        valid = ~isnan(ap(i,:)); % Find valid (non-NaN) data points
        if sum(valid) > 1  % Ensure at least two valid points
            yiyi(i, :) = spline(wl(valid), ap(i, valid), xixi);
        end
    end

    % Assign to absspec
    absspec = zeros(size(ap,1), length(wavelength));
    absspec(:, minwavel*10:maxwavel*10) = yiyi; % Now correctly sized

    absspec(:, 1:minwavel*10-1) = ones(1, size(1:minwavel*10-1,2)) .* absspec(:, minwavel*10);
    aspecprime = absspec';

    meassignal6 = NaN(size(aspecprime, 2), size(wl, 2));
    for i = 1:size(aspecprime, 2)
        measur2 = aspecprime(:,i) .* filtfunc; % the measured signal for every filter factor.
        meassignal6(i,:) = 0.1 * sum(measur2); % The measured spectrum at a wavelength i is the sum of what a filter measured at
    end
    abscorr = ap - meassignal6 + ap;
    fprintf('Done\n')
else
    fprintf('Skip absorption spectra correction\n')
    abscorr = aph;

end

% Peak center values ("peak_loc") determined using a interative
% approach that allows the location to vary (uses the matlab
% function LSQNONLIN), and are rounded to nearest integer.
% Sigma values ("lsqsig") are determined similarly. FWHM = sigma*2.355
peak_loc=[434,453,470,492,523,550,584,617,638,660,675];

lsqsig=[12,12,13,16,14,14,16,13,11,11,10];

onenm = 410:1:700;

fprintf('Gaussian decomposition ... ')

% Preallocate array to store interpolated spectra
acorr2onenm = NaN(size(abscorr, 1), numel(onenm));

% Loop over each row of abscorr
for i = 1:size(abscorr, 1)
    % Find non-NaN indices for the current row
    non_nan_indices = ~isnan(abscorr(i, :));
    
    % Extract non-NaN wavelengths and corresponding absorption values
    wavelengths = wl(non_nan_indices);
    absorption_values = abscorr(i, non_nan_indices);
    
    % Check if absorption_values is not empty
    if ~isempty(absorption_values) && length(absorption_values) > 1
        % Interpolate the absorption spectrum to onenm
        interpolated_spectrum = interp1(wavelengths, absorption_values, onenm, 'spline');

        % Store the interpolated spectrum
        acorr2onenm(i, :) = interpolated_spectrum;
    else
        acorr2onenm(i, :) = NaN;
    end

    acorr2onenm(acorr2onenm(i, :) < 0) = NaN;
end


%define the matrix of component Gaussian functions using the peaks
%and widths (sigma) above
coef2 = exp(-0.5 .* (((onenm .* ones(size(peak_loc,2),1))' - peak_loc .* ...
    ones(size(onenm,2),1)) ./ lsqsig) .^ 2);

% Preallocate array to store interpolated uncertainties
apunc_int = NaN(size(apunc, 1), numel(onenm));

% Loop over each row of apunc
for i = 1:size(apunc, 1)
    % Find non-NaN indices for the current row
    non_nan_indices = ~isnan(apunc(i, :));

    % Extract non-NaN wavelengths and corresponding uncertainty values
    wavelengths = wlunc(non_nan_indices);
    uncertainty_values = apunc(i, non_nan_indices);

    % Check if uncertainty_values is not empty
    if ~isempty(uncertainty_values) && length(uncertainty_values) > 1
        % Interpolate the uncertainty spectrum to onenm
        interpolated_uncertainty = interp1(wavelengths, uncertainty_values, onenm, 'spline');

        % Store the interpolated uncertainty
        apunc_int(i, :) = interpolated_uncertainty;
    else
        apunc_int(i, :) = NaN;
    end

    apunc_int(apunc_int(i, :) < 0) = NaN;
end

% Compute the mean uncertainty spectrum across all stations (ignoring NaNs)
mean_apunc = nanmean(apunc_int, 1);  % [1 × 321]

% Now apunc_int contains the interpolated uncertainties for all stations
amps = NaN(size(aph,1), size(coef2,2));
sumspec_temp = NaN(size(coef2,1), size(aph,1));
compspec_temp = NaN(size(coef2,1), size(coef2,2), size(aph,1));

for i = 1:size(acorr2onenm, 1)

    % If ap exists but uncertainty is all NaN, replace with mean
    if all(isnan(apunc_int(i, :))) && any(~isnan(acorr2onenm(i, :)))
        apunc_int(i, :) = mean_apunc;
    end

    % Identify valid wavelength indices (non-NaN in both signal and uncertainty)
    valid_idx = ~isnan(acorr2onenm(i, :)) & ~isnan(apunc_int(i, :));
    
    if sum(valid_idx) < size(coef2, 2)
        % Not enough valid data points to resolve all 13 components
        continue
    end

    % Get valid ap spectrum and uncertainties
    target = acorr2onenm(i, valid_idx);            % [n_valid × 1]
    unc = apunc_int(i, valid_idx);                 % [n_valid × 1]
    
    % Weight both the model and the target by uncertainty
    target_w = target ./ unc;                      % [n_valid × 1]
    coef2_w = coef2(valid_idx, :) ./ unc';         % [n_valid × 13]

    % Solve the non-negative least squares problem
    amps(i, :) = lsqnonneg(coef2_w, target_w');    % [1 × 13]

    if all(amps(i, :) == 0)
        amps(i, :) = NaN;  % No meaningful fit
        continue
    end

    % Reconstruct sum and component spectra
    sumspec_temp(:, i) = sum(amps(i, :) .* coef2, 2);               % [321 × 1]
    compspec_temp(:, :, i) = amps(i, :) .* coef2;                   % [321 × 13]
end

valid_station_mask = ~all(isnan(amps), 2);
valid_indices = find(valid_station_mask);

n_wl = length(wl);
n_comp = size(coef2,2);
n_stations = size(aph,1);

compspec = NaN(n_wl, n_comp, n_stations);
sumspec = NaN(n_stations, n_wl);

for idx = 1:length(valid_indices)
    i = valid_indices(idx);

    compspec(:, :, i) = interp1(onenm', squeeze(compspec_temp(:,:,i)), wl, 'spline');
    sumspec(i, :) = interp1(onenm, sumspec_temp(:,i), wl, 'spline');
end


fprintf('Done\n')

end
