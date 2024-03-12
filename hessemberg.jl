begin
	
function isHeisenberg(H)
	n, m = size(H);

	if (n != m)
		return false;
	end

	isUpper = true;
	isLower = true;

	for i = 1:n
		for j = (i + 2):m
			if (H[i, j] != 0.0)
				isUpper = false;
			end

			if (H[j, i] != 0.0)
				isLower = false;
			end
		end

		if (!(isUpper || isLower))
			break;
		end
		
	end

	return isUpper || isLower;
	
end

# testing
H = [6.7 2.3 0.0 0.0 1.0;
	 9.1 4.5 4.0 0.0 0.0;
	 0.0 7.2 3.0 8.3 0.0;
	 0.0 1.0 0.3 1.6 1.1;
	 0.0 0.0 0.0 2.9 0.0];

println(isHeisenberg(H)); # false
		
end
