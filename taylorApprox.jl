begin

using Plots;

function F(x::Real, n::Integer)
	f(x) = exp(2x) + exp(4x / 3);
	fDerivative(x, p) = 2^p * exp(2x) + (4/3)^p * exp(4/3 * x);

	y = f(0);

	# instead of storing x^i and i!, wich would lead to too large numbers,
	# store their division directly, leading to more numerical stability
	xPowerOverFactorial = 1; # x^0 / 0!

	for i = 1:n
		xPowerOverFactorial *= x / i;
		y += fDerivative(0, i) * xPowerOverFactorial;
	end

	return y;
end


function numeroTermos(x::Real, atol::Real)
	f(x) = exp(2x) + exp(4x / 3);

	l = 0;
	r = 1;
	while (!isapprox(F(x, r), f(x), atol = atol))
		l = r
		r *= 2;
	end

	while (l < r - 1)
		m = trunc(Int32, l + (r - l) * 0.5);

		if (!isapprox(F(x, m), f(x), atol = atol))
			l = m;
		else
			r = m
		end
	end
	
	return r;
end




plot(x -> numeroTermos(x, 1e-3), -5, 5);
  

end
