% This function calculates oscillation index (oi) using Gobbel and Fike's formula.
% Ref: Reference [19] of Wu, 2003 doi/10.1002/mrm.10522/
% Input: residue vector scaled by CBF
% Output: oscillation index value

function oscillation_index = calculate_oi_Gobbel_Fike(residue_vector)
	sum_absolute = 0; % Initiate value to save summation of absolute values

	for k = 3 : length(residue_vector)
		sum_absolute = sum_absolute + abs(residue_vector(k) - 2 * residue_vector(k - 1) + residue_vector(k - 2));
	end

	oscillation_index = (1 / length(residue_vector)) * (1 / max(residue_vector)) * sum_absolute;

end
