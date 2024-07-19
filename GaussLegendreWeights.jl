begin
	using LinearAlgebra;
	using Plots;

  
  function multPolynomials(coeffs1, coeffs2)
  	N = length(coeffs1);
  	M = length(coeffs2);
  
  	coeffs = zeros(N + M - 1);
  
  	for i = 1:N
  		for j = 1:M
  			coeffs[i + j - 1] += coeffs1[i] * coeffs2[j];
  		end
  	end
  
  	return coeffs;
  end

  # Horner's method for evaluating a polynomial
  function evaluatePolynomial(coeffs, x)
  	v = coeffs[end];
  
  	for i = length(coeffs)-1:-1:1
  		v = v * x + coeffs[i];
  	end
  
  	return v;
  end

  # evaluate the derivative in-place
  function evaluatePolynomialDerivative(coeffs, x)
  	v = coeffs[end] * (length(coeffs) - 1);
  
  	for i = length(coeffs)-1:-1:2
  		v = v * x + coeffs[i] * (i - 1);
  	end
  
  	return v;
  end

  # the eigenvalues of a polynomial's companion matrix are it's roots
  function polynomialRoots(coeffs)
  	N = length(coeffs) - 1;
  	C = zeros(N, N);
  
  	C[1, end] = -coeffs[1] / coeffs[end];
  	for i = 2:N
  		C[i, end] = -coeffs[i] / coeffs[end];
  		C[i, i - 1] = 1.0;
  	end
  
  	return eigvals(C);
  end

  # returns the coefficients of the n-th Legendre polynomial
  function LegendreCoefficients(n)
  	coeffs = [1];
  	for i = 1:n
  		coeffs = multPolynomials(coeffs, [-1, 0, 1]); # multiply by -1 + x²
  	end

    # derive n times
  	for i = 1:n
  		for j = 2:length(coeffs)+1-i
  			coeffs[j - 1] = coeffs[j] * (j - 1);
  		end
  	end
  
  	resize!(coeffs2, n + 1);
  
  	for j = 1:length(coeffs)
  		coeffs[j] /= (2.0 ^ n) * factorial(n);
  	end
  
  	return coeffs;
  	
  end


  # returns all the points and their correspondant weights for gaussian quadrature
  function GaussQuadraturePointsAndWeights(n)
  
  	points = [];
  	weights = [];
  
  	coeffs = LegendreCoefficients(n);
  	r = polynomialRoots(coeffs);
  
  	for i = 1:length(r)
  		x = real(r[i]);
  		xd = evaluatePolynomialDerivative(coeffs, x);
  		
  		push!(points, x)
  		push!(weights, (2 / ((1 - x * x) * xd * xd)))
  	end
  
  	return points, weights;
  end

  
end