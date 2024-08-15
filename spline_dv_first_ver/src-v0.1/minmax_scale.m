function [xyz_scal, xyzp_scal] = minmax_scale(xyz, xyzp)
    % Re-scaling mean and CV, keeping dropout rate intact
    % xyz component 1 must be mean
    % xyz component 2 must be CV
    % xyz component 3 must be dropuot rate

    max_x = max( xyz(:,1) );
    max_y = max( xyz(:,2) );

    xyz_scal = [xyz(:,1)./max_x xyz(:,2)./max_y xyz(:,3) ];
    xyzp_scal = [xyzp(:,1)./max_x xyzp(:,2)./max_y xyzp(:,3)];
end