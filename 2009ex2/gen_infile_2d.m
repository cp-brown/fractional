function gen_infile_2d(name, method, difcoef, order, focus, t1, t2, nz, theta, dx, dy, dt, ax, bx, nx, ay, by, ny, t)

    fid = fopen(name, 'w');
    fprintf(fid, '%d', method);
    fprintf(fid, '    method\n');
    fprintf(fid, '%.15g', difcoef);
    fprintf(fid, '    diffusion coefficient\n');
    fprintf(fid, '%.15g', order);
    fprintf(fid, '    order\n');
    fprintf(fid, '%.15g', focus);
    fprintf(fid, '    focus\n');
    fprintf(fid, '%.15g', t1);
    fprintf(fid, '    t1\n');
    fprintf(fid, '%.15g', t2);
    fprintf(fid, '    t2\n');
    fprintf(fid, '%d', nz);
    fprintf(fid, '    nz\n');
    fprintf(fid, '%.15g', theta);
    fprintf(fid, '    theta\n');
    fprintf(fid, '%.15g', dx);
    fprintf(fid, '    dx\n');
    fprintf(fid, '%.15g', dy);
    fprintf(fid, '    dy\n');
    fprintf(fid, '%.15g', dt);
    fprintf(fid, '    dt\n');
    fprintf(fid, '%d', nx);
    fprintf(fid, '    nx\n');
    fprintf(fid, '%d', ny);
    fprintf(fid, '    ny\n');
    fprintf(fid, '%d', length(t));
    fprintf(fid, '    nt\n');
    fprintf(fid, '%.15g ', linspace(ax, bx, nx));
    fprintf(fid, '\n');
    fprintf(fid, '%.15g ', linspace(ay, by, ny));
    fprintf(fid, '\n');
    fprintf(fid, '%.15g ', t);
    fclose(fid);

end