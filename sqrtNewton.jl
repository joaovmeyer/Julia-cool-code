function sqrtNewton(x::Real, n::Integer)
	y = x * 0.5;
	
	for i = 1:n
		# add eps(y) for numerical stability
		y = (y + x / (y + eps(y))) * 0.5;
	end

	return y;
end
