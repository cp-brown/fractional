exact = ex1_2009_exact(2,1);

% Table 1
N = [10, 20, 30, 40, 60, 80, 100, 120];
err = NaN(length(N), 2);
runtime = 0;
for j = 1 : length(N)

    file = ['out/m1o0nz', num2str(N(j)), '.txt'];
    fid = fopen(file);
    runtime = runtime + str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    err(j,1) = abs(u - exact);
    fclose(fid);

    file = ['out/m1o2nz', num2str(N(j)), '.txt'];
    fid = fopen(file);
    runtime = runtime + str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    err(j,2) = abs(u - exact);
    fclose(fid);

end
disp('2009 Paper, Table 1:');
disp('  N      err 0      err 2');
disp('-------------------------');
for j = 1 : length(N)
    fprintf('%3d   %.2e   %.2e\n', N(j), err(j,1), err(j,2));
end
fprintf('Total computation time: %.3e\n', runtime);
disp(' ');


% Table 2
nz = [10, 20, 40, 60, 80, 100, 120, 160, 200];
err = NaN(length(nz), 2);
runtime = 0;
for j = 1 : length(nz)

    file = ['out/m2nz', num2str(nz(j)), '.txt'];
    fid = fopen(file);
    runtime = runtime + str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    err(j,1) = abs(u - exact);
    fclose(fid);

    file = ['out/m3nz', num2str(nz(j)), '.txt'];
    fid = fopen(file);
    runtime = runtime + str2num(fgets(fid));
    t = str2num(fgets(fid));
    u = str2num(fgets(fid));
    err(j,2) = abs(u - exact);
    fclose(fid);

end
disp('2009 Paper, Table 2:');
disp('  N     err U0     err U~');
disp('-------------------------');
for j = 1 : length(nz)
    fprintf('%3d   %.2e   %.2e\n', nz(j), err(j,1), err(j,2));
end
fprintf('Total computation time: %.3e\n', runtime);