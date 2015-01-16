% This function performs deconvolution using block circulant matrix
% Input parameters:
% signal_vector: ASL signal time series of one voxel
% aif_vector: arterial input function of the corresponding ASl signal time series
% deltaTI: delta_ti
% Ouput: residue vector (also time series) scaled by CBF
% Ref: Deconvolution Using a Block-Circulant Matrix, Wu 2003, doi/10.1002/mrm.10522

function residue_scaled_vector = svd_block_circulant(signal_vector, aif_vector, deltaTI)

	n_ti = size(signal_vector); % get the dimension of signal_matrix, for now n_column = n_ti

	% Perform zero padding
	% By zero-padding the N-point time series Ca(t) and C(t) to length L, where L ≥ 2N,
	% time aliasing can be avoided
	
	aif_scaled_vector = deltaTI * aif_vector; % Scale the AIF vector by deltaTI
	% Create a lower triangular matrix in the form of matrix A (Wu, 2003)
	aif_triangular_matrix = convert_to_low_tri(aif_scaled_vector);


end

