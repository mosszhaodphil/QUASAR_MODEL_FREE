% This function performs deconvolution using block circulant matrix
% Input parameters:
% signal_vector: ASL signal time series of one voxel
% aif_vector: arterial input function of the corresponding ASl signal time series
% deltaTI: delta_ti
% Output: residue vector (also time series) scaled by CBF
% Ref: Deconvolution Using a Block-Circulant Matrix, Wu 2003, doi/10.1002/mrm.10522

function residue_scaled_vector = svd_block_circulant(signal_vector, aif_vector, deltaTI)

	n_ti = length(signal_vector); % get the dimension of signal_matrix, for now n_column = n_ti

	% Perform zero padding
	% By zero-padding the N-point time series Ca(t) and C(t) to length L, where L ≥ 2N,
	% time aliasing can be avoided
	new_length = floor(n_ti * 2.2); % now new_length ≥ 2 * n_ti

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
	[U, S, V_transpose] = svd(aif_block_circulant_matrix);

	% Construct W, V, U_transpose such that V * W * U_transpose = D or inverse of A (Wu, 2003)
	W           = inv(S);
	V           = (V_transpose);
	U_transpose = (U);

	% According to (Wu, 2003), we need a threshold p_SVD such that
	% if value of S is less than p_SVD, we set the corresponding value in W to zero.
	% Since S is decreasing along the diagonal,
	% we can set p_SVD to be the smallest value along the diagonal of S.
	% As such, p_SVD will increase along the diagonal of S (from bottom to top).
	% As a result, the values along the diagonal of W will become zero from bottom to top.

	% We begin by using all singular values (diagonal) of W
	residue_scaled_vector = V * W * U_transpose * signal_padding_vector; % calculate residue scaled by CBF
	oi = calculate_oi_Gobbel_Fike(residue_scaled_vector); % Calculate the current oscillation_index(oi)
	oi_threshold = 0.1; % set an oi threshold to be updated

	% Start removing singular values of W one by one from bottom to top
	% Do so until the calculated oi is less than the threshold
	j = length(residue_scaled_vector);
	while((oi > oi_threshold) && j > 1)
		p_SVD = W(j, j);
		% Set the lowest non-zero singular value of W to zero
		W(j, j) = 0; % This W(j, j) corresponds to the smallest S(j, j). We set it to zero.

		% Update residue vector and oi
		residue_scaled_vector = V * W * U_transpose * signal_padding_vector; % calculate residue scaled by CBF
		oi = calculate_oi_Gobbel_Fike(residue_scaled_vector); % Calculate the current oscillation_index(oi)
		j = j - 1;
	end

	% Now the residue vector (scaled by CBF) should be the output
	% This length of residue vector is new_length or (n_ti + padding)
	% CBF value is the largest element of this residue vector

end

