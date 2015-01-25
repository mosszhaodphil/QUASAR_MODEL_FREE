% This function correct the calculated CBF (SVD WU) to absolute scale
% Input:
% unscaled_cbf_value: unscaled CBF value, the maximum element of residue vector
% Output:
% corrected_cbf_value: corrected CBF value

function corrected_cbf_value = correct_cbf_scale(unscaled_cbf_value)

	corrected_cbf_value = unscaled_cbf_value * 6000 / 2 / 0.91;

end

