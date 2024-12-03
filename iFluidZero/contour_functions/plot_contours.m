function plot_contours(contour, contour_color, fill_color)

    hold on        
    
    for i = 1:size(contour, 3)
        Gi = contour(:,:,i);
        Gi(end+1,:) = Gi(1,:); % close the loop
        
        if nargin == 1
            plot(Gi(:,1), Gi(:,2), '.-')
        elseif nargin == 2
            plot(Gi(:,1), Gi(:,2), '.-', 'color', contour_color)
        else
            fill(Gi(:,1),Gi(:,2),fill_color)
            plot(Gi(:,1), Gi(:,2), '.-', 'color', contour_color)
        end
        
    end
    
    hold off
end