% This function convert the input matrix into a block circulant matrix (actually just circular matrix)
% Ref: Deconvolution Using a Block-Circulant Matrix (D) of 'Wu 2003, doi/10.1002/mrm.10522'

function block_circulant_matrix = create_block_circulant_matrix(input_matrix)

	[x, y] = size(input_matrix); % input matrix should be square so x = y
	block_circulant_matrix = zeros(x, y); % create an empty matrix of the same size with input_matrix

	% Construct block circulant matrix
	% replacing matrix A with a block-circulant matrix D,
	% whose elements are d(i,j) = a(i,j) for j ≤ i,
	% and d(i,j) = a(L+i-j,0) otherwise. L is size of matrix a
	for i = 1 : x
		for j = 1 : y
			if((j < i) || (j == i))
				block_circulant_matrix(i, j) = input_matrix(i, j);
			else
				block_circulant_matrix(i, j) = input_matrix((x + 1) + i - j, 1); % The index of input_matrix starts from 1 not zero, so we need plus one
		end
	end

end

