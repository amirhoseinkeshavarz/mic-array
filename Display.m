%% plot
faceopaque = false;
icahedron_plot = true;

if icahedron_plot == true
    x = spacePointsFine(:, 1);
    y = spacePointsFine(:, 2);
    z = spacePointsFine(:, 3);
    if faceopaque == true
        h = trisurf(fMatFine,x,y,z, 'FaceColor', 'white', 'EdgeColor', 0.2*[1 1 1], 'LineWidth', 1, 'FaceAlpha', 1);
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black', 'DisplayName', 'discrete grid');
        plot3(micPosition(:, 1), micPosition(:, 2), micPosition(:, 3), 'sg', 'DisplayName', 'microphones')
    else
        h=trisurf(fMatFine,x,y,z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1],'LineWidth', 1 );
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black', 'DisplayName', 'discrete grid');
        plot3(micPosition(:, 1), micPosition(:, 2), micPosition(:, 3), 'sg', 'MarkerFaceColor', 'green', 'DisplayName', 'microphones')
    end
    
    xlabel('x-axis (m)');
    ylabel('y-axis (m)');
    zlabel('z-axis (m)');
    title( ['test points']);
    
    axis equal;
    grid off;
    
end