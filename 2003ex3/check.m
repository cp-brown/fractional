exact = @(x,t) sin(x) * (1+t) .* exp(-t) + sin(2*x) * cos(t) .* exp(-2*t);

% Spatial convergence
nt = 12;
nx = [10, 20, 40, 80];
err = NaN(nt, length(nx));
for j = 1 : length(nx)

    file = ['out/sp', num2str(nx(j)), '.txt'];
    fid = fopen(file);
    x = str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    fclose(fid);
    u = reshape(u, length(x), length(t));
    u = [zeros(1,length(t)); u; zeros(1,length(t))];
    x = x(:);
    x = [0; x; pi];

    for k = 1 : length(t)
        err(k,j) = integral(@(s)abs(interp1(x, u(:,k), s) - exact(s,t(k))), ...
        x(1), x(end));
    end

end
roc = log2( err(:,1:end-1) ./ err(:,2:end) );
disp('2003 Paper, Table 8, method 1:')
disp('  t     err_10     err_20    r_20     err_40    r_40     err_80    r_80');
disp('---------------------------------------------------------------------------');
for j = 1 : length(t)
    fprintf('%.1f   %.2e   %.2e   %5.2f   %.2e   %5.2f   %.2e   %5.2f\n', ...
        t(j), err(j,1), err(j,2), roc(j,1), err(j,3), roc(j,2), ...
        err(j,4), roc(j,3));
end
disp(' ');


% Temporal convergence
nz = [10, 20, 30, 40, 60, 80, 100, 120];
err = NaN(length(nz), 3);
for j = 1 : length(nz)

    file = ['out/tm1n', num2str(nz(j)), '.txt'];
    fid = fopen(file);
    x = str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    fclose(fid);
    u = u(:);
    u = [0; u; 0];
    x = x(:);
    x = [0; x; pi];
    err(j,1) = integral(@(s)abs(interp1(x, u, s) - exact(s,t)), x(1), x(end));

    file = ['out/tm2n', num2str(nz(j)), '.txt'];
    fid = fopen(file);
    x = str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    fclose(fid);
    u = u(:);
    u = [0; u; 0];
    x = x(:);
    x = [0; x; pi];
    err(j,2) = integral(@(s)abs(interp1(x, u, s) - exact(s,t)), x(1), x(end));

    file = ['out/tm3n', num2str(nz(j)), '.txt'];
    fid = fopen(file);
    x = str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    fclose(fid);
    u = u(:);
    u = [0; u; 0];
    x = x(:);
    x = [0; x; pi];
    err(j,3) = integral(@(s)abs(interp1(x, u, s) - exact(s,t)), x(1), x(end));
    
end
disp('Temporal convergence:');
disp('  N      err 1      err 2      err 3');
disp('------------------------------------');
for j = 1 : length(nz)
    fprintf('%3d   %.2e   %.2e   %.2e\n', nz(j), err(j,1), err(j,2), err(j,3));
end