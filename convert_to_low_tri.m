% This code is the implementation the following papers
% MA Chappell (2012) doi: 10.1002/mrm.24372 (MACQ)
% MA Chappell (2012) doi: 10.1002/mrm.24260 (MACD)
% ET Petersen (2006) doi: 10.1002/mrm.20784 (ETP)
% M  Günther  (1998) doi: 10.1002/mrm.1284 (MG)

% This function converts a vector to 'lower triangluar matrix' as AIF in eq [12] (ETP)

function lower_tri_matrix = convert_to_low_tri(v)

	lower_tri_matrix = zeros(length(v), length(v)); % create an empty square matrix

	for k = 1 : length(v)
		for j = 1 : length(v)
			if(j >= k)
				lower_tri_matrix(j, k) = v(j - k + 1);
			else
				% lower_tri_matrix(j, k) remains zero
			end
		end
	end
end