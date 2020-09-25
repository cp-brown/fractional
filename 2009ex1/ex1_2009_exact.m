function u = exact(t, a)
% Exact solution for 2009 example 1

    E = @(t) exp(t.^2) .* erfc(-t);

    f = @(x) exp(-x).*cos(pi()*x);
    u = NaN(size(t));
    for j = 1:length(t)
        u(j) = E(-a*sqrt(t(j))) + integral(@(y) E(-a*sqrt(t(j))*y) ...
            .* f(t(j)-t(j)*y.^2) * 2 * t(j) .* y, 0, 1);
    end
end
