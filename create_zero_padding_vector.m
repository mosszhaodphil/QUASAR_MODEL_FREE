% This function creates a zero padding vector
% Zeros will be padded at the end of the input vector
% Input parameters:
% input_vector: vector to be padded
% new_length: length of the zero padding vector
% Output: padding vector
function padding_vector = create_zero_padding_vector(input_vector, new_length)

	padding_vector = zeros(new_length, 1); % Create an empty vector of size new_length

	% copy the first n element of input_vector to padding vector
	% n is the size of input_vector
	for j = 1 : length(input_vector)
		padding_vector(j) = input_vector(j);
	end

end
