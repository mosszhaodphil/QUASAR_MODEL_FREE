% This function performs model free analysis to estimate CBF with block circulant SVD
% Input:
% tissue_matrix: 4D QUASAR (tissue) signal (Control - Label)
% aif_matrix: 4D arterial input function (aif) signal
% abv_matrix: 3D arterial blood volume
% deltaTI: sampling rate (default: 0.3s)
% Output:
% cbf_matrix: 3D estimated CBF values
% residue_matrix: 4D residue result from block circulant SVD
% Reference: Deconvolution Using a Block-Circulant Matrix, Wu 2003, doi/10.1002/mrm.10522

function [cbf_matrix residue_matrix] = calculate_cbf_model_free(tissue_matrix, aif_matrix, abv_matrix, deltaTI)

	[x, y, z, t] = size(tissue_matrix); % get the dimension of input matrix

	%residue_matrix = zeros(x, y, z, t);
	cbf_matrix     = zeros(x, y, z);

	for k = 1 : z
		for j = 1 : y
			for i = 1 : x
				% Get current aif (time series) of current voxel
				aif_vector = aif_matrix(i, j, k, :);
				aif_vector = aif_vector(:);

				% Get current tissue signal (time series) of current voxel
				tissue_vector = tissue_matrix(i, j, k, :);
				tissue_vector = tissue_vector(:);

				% Now use block circulant SVD to perform model free analysis
				residue_vector      = svd_block_circulant(tissue_matrix(i, j, k, :), aif_matrix(i, j, k, :), deltaTI); % Get residue vector
				unscaled_cbf_value  = calculate_perfusion_from_residue_vector(residue_vector); % Get CBF value of current residue vector
				corrected_cbf_value = correct_cbf_scale(unscaled_cbf_value); % Correct CBF into absolute scale

				% Store the result into matrix
				% Note: dimension of residue_matrix is (x, y, z, t + padding)
				% padding is the number of padded zeros added in block circulant SVD
				cbf_matrix(i, j, k) = corrected_cbf_value;
				residue_matrix(i, j, k, :) = residue_vector;
			end
		end
	end


end

