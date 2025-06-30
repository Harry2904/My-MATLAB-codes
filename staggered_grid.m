function staggered_grid(imax, jmax)
    % Domain size (Lx = Ly = 1 for non-dimensional form)
    Lx = 1; 
    Ly = 1;
    dx = Lx / (imax - 2);
    dy = Ly / (jmax - 2);
    
    %grid Points
    [Xp, Yp] = meshgrid(linspace(dx/2, Lx-dx/2, imax), linspace(dy/2, Ly-dy/2, jmax));
    [Xu, Yu] = meshgrid(linspace(0, Lx, imax+1), linspace(dy/2, Ly-dy/2, jmax));
    [Xv, Yv] = meshgrid(linspace(dx/2, Lx-dx/2, imax), linspace(0, Ly, jmax+1));
    
    % Plot staggered grid
    figure;
    hold on;
    plot(Xp, Yp, 'ro', 'MarkerSize', 6, 'DisplayName', 'Pressure (P)');
    plot(Xu, Yu, 'bs', 'MarkerSize', 6, 'DisplayName', 'U-velocity (U)');
    plot(Xv, Yv, 'g^', 'MarkerSize', 6, 'DisplayName', 'V-velocity (V)');
    
    % Grid lines
    for i = 1:imax
        xline(Xp(1, i), '--k');
    end
    for j = 1:jmax
        yline(Yp(j, 1), '--k');
    end
    
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Staggered Grid Representation');
    axis equal;
    grid on;
    hold off;
end

