nz = [10, 20, 30, 40, 60, 80, 100, 120];
for j = 1 : length(nz)

    file = ['in/tm1n', num2str(nz(j)), '.txt'];
    gen_infile_2d(file, 1, 1, -0.5, 1, 0.5, 5, nz(j), 0.1, 0.001, 0.001, 0.001, 0, 4, 50, 0, 4, 50, 2);

    file = ['in/tm2n', num2str(nz(j)), '.txt'];
    gen_infile_2d(file, 2, 1, -0.5, 1, 0.5, 5, nz(j), 1, 0.001, 0.001, 0.001, 0, 4, 50, 0, 4, 50, 2);

    file = ['in/tm3n', num2str(nz(j)), '.txt'];
    gen_infile_2d(file, 3, 1, -0.5, 0, 0.5, 5, nz(j), 1, 0.001, 0.001, 0.001, 0, 4, 50, 0, 4, 50, 2);

end
