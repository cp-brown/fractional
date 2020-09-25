function u = exact(x, y, t, homogeneous)
% Exact solution for 2009 example 2

    E = @(t) exp(t.^2) .* erfc(-t);

    phi = @(x,y,j,k) sin(j*pi*x(:)/4) * sin(k*pi*y(:)'/4);

    % Compute u11
    if not(homogeneous)
        f = @(s) exp(-s/2).*cos(s);
    else
        f = @(s) zeros(size(s));
    end
    a = pi^2/16 * (1+1);
    u11 = NaN(size(t));
    for j = 1:length(t)
        u11(j) = E(-a*sqrt(t(j))) + integral(@(z) E(-a*sqrt(t(j))*z) .* f(t(j)-t(j)*z.^2) * 2 * t(j) .* z, 0, 1);
    end

    % Compute u21
    if not(homogeneous)
        f = @(s) 0.5*exp(-s).*cos(pi*s);
    else
        f = @(s) zeros(size(s));
    end
    a = pi^2/16 * (4+1);
    u21 = NaN(size(t));
    for j = 1:length(t)
        u21(j) = E(-a*sqrt(t(j))) + integral(@(z) E(-a*sqrt(t(j))*z) .* f(t(j)-t(j)*z.^2) * 2 * t(j) .* z, 0, 1);
    end

    u = NaN(length(x), length(y), length(t));
    for j = 1 : length(t)
        u(:,:,j) = u11(j) * phi(x,y,1,1) - u21(j) * phi(x,y,2,1);
    end
end
