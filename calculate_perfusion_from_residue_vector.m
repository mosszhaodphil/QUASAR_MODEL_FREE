% This function calculate the perfusion value from residue vector
% Input parameters:
% residue_vector: residue vector result from svd_block_circulant deconvolution by Wu
% perfusion: perfusion value
% Ref: Perfusion is the maximum value of residue vector, Wu 2003, doi/10.1002/mrm.10522

function perfusion = calculate_perfusion_from_residue_vector(residue_vector)

	perfusion = max(residue_vector);
	
end

