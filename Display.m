%% plot

if icahedron_plot == true
    x = vMat(:, 1);
    y = vMat(:, 2);
    z = vMat(:, 3);
    if faceopaque == true
        h=trisurf(fMat,x,y,z, 'FaceColor', 'white', 'EdgeColor', 0.2*[1 1 1], 'LineWidth', 1, 'FaceAlpha', 1);
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black');
        plot3(mics_locs(:, 1), mics_locs(:, 2), mics_locs(:, 3), 'sg')
    else
        h=trisurf(fMat,x,y,z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1],'LineWidth', 1 );
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black');
        plot3(mics_locs(:, 1), mics_locs(:, 2), mics_locs(:, 3), 'sg', 'MarkerFaceColor', 'green')
    end
    
    xlabel('x-axis (m)');
    ylabel('y-axis (m)');
    zlabel('z-axis (m)');
    title( ['test points']);
    
    axis equal;
    grid off;
    
end