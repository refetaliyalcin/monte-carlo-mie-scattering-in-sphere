function [s_x_,s_y_,s_z_]=scatter_mc(cos_theta,s_x,s_y,s_z)

    sin_theta=sqrt(1 - cos_theta*cos_theta);
    
    R_azimuth=rand();
    phi=2*pi*R_azimuth;
    cos_phi=cos(phi);
    if phi<pi
        sin_phi=sqrt(1-cos_phi*cos_phi);
    else
        sin_phi=-sqrt(1-cos_phi*cos_phi);
    end
    
    if (s_z==1)
        s_x_ = sin_theta*cos_phi;
        s_y_ = sin_theta*sin_phi;
        s_z_ = cos_theta;                                           
    elseif (s_z==-1)
        s_x_ = sin_theta*cos_phi;
        s_y_ = sin_theta*sin_phi;
        s_z_ = -cos_theta;
    else     
        denom = sqrt(1 - s_z*s_z);
        s_x_ = sin_theta*(s_x * s_z * cos_phi - s_y * sin_phi) / denom + s_x * cos_theta;
        s_y_ = sin_theta*(s_y * s_z * cos_phi + s_x * sin_phi) / denom + s_y * cos_theta;
        s_z_ = -denom*sin_theta*cos_phi + s_z*cos_theta;
    end
end