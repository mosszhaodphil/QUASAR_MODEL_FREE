% This function performs deconvolution using block circulant matrix
% Input parameters:
% signal_vector: ASL signal time series of one voxel
% aif_vector: arterial input function of the corresponding ASl signal time series
% deltaTI: delta_ti
% Ouput: residue vector (also time series) scaled by CBF
% Ref: Deconvolution Using a Block-Circulant Matrix, Wu 2003, doi/10.1002/mrm.10522

function residue_scaled_vector = svd_block_circulant(signal_vector, aif_vector, deltaTI)

	n_ti = length(signal_vector); % get the dimension of signal_matrix, for now n_column = n_ti

	% Perform zero padding
	% By zero-padding the N-point time series Ca(t) and C(t) to length L, where L ≥ 2N,
	% time aliasing can be avoided
	new_length = floor(n_ti * 2.1); % now new_length ≥ 2 * n_ti

	% Create zero padding vectors for signal and aif
	aif_padding_vector    = create_zero_padding_vector(aif_vector, new_length);
	signal_padding_vector = create_zero_padding_vector(signal_vector, new_length);
	
	% Scale the AIF padding vector by deltaTI
	aif_scaled_vector = deltaTI * aif_padding_vector; 

	% Create a lower triangular matrix in the form of matrix A (Wu, 2003)
	aif_triangular_matrix = convert_to_low_tri(aif_scaled_vector); % aif_triangular_matrix is the matrix A in (Wu, 2003)

	% Create block circulant matrix (D in Wu 2003)
	aif_block_circulant_matrix = create_block_circulant_matrix(aif_triangular_matrix);

	% Perform singular value decomposition of the inverse of aif_block_circulant_matrix
	[V, W, U_transpose] = svd(inv(aif_block_circulant_matrix));

	residue_scaled_vector = V * W * U_transpose * aif_padding_vector;

	% Now perform oscillation index by Gobbel and Fike

end

