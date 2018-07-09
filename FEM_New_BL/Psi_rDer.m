function val = Psi_rDer(r, t, epsilon)
%PSI_RDER The r derivative of the Psi function

val = -1 / (4 * pi^(1/2) * epsilon^(3/2)) .* ...
    exp(-r.^2 ./ (4 .* epsilon .* t)) .* ...
    ( r.^3 ./ (4 .* epsilon.^2 .* t.^(7/2)) - (3/2) .* r ./ (epsilon .* t.^(5/2)) );

end

