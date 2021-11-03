h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for time = Dt:10*Dt:TOTAL_TIME

    T_plot = TMID(:,round(time/Dt));
    plot(x',T_plot); title(sprintf("Time = %10.3e", time));
    xlim([0, max(x)]); grid on
%     ylim([421.5, 424.5]);
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
%     if time == TOTAL_TIME
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
end