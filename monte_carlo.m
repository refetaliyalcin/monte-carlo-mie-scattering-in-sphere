function [r_tot,t_tot,a_tot,t_direct_tot,acilar,donme_sayisi]=monte_carlo(photon_number,R,scat_prob,mu_tot,inv_cdf_cos,n_cdf_random,g,use_HG)
r_tot=0;
t_tot=0;
t_direct_tot=0;
a_tot=0;
donme_tot=0;
acilar=zeros(180,1);
for i=1:photon_number
    donme_no=0;
    r_no=0;
    t_no=0;
    t_direct_no=0;
    a_no=0;
    ever_extinct=0;
    alive=1;
    outside=1;
    while outside
        x=R*(1-2*rand()); %x component of position vector
        y=R*(1-2*rand()); %y component of position vector
        if sqrt(x*x+y*y) <= R
             outside=0;
        end
    end
    z=-sqrt(R*R-x*x-y*y);%z component of position vector

    s_x=0;%x component of direction vector
    s_y=0;%y component of direction vector. could be defined better along with s_x but it is a 1D code anyway, nothing will change
    s_z=1; %z component of direction vector 
    l_beta=-log(rand())/mu_tot; %excitation length
    while alive   
        intersections=intersectLineSphere([x y z s_x s_y s_z],[0 0 0 R]);
        point_1=intersections(1,:);
        point_2=intersections(2,:);
        if point_1(3)>point_2(3)
            point_up=point_1;
            point_down=point_2;
        else
            point_up=point_2;
            point_down=point_1;
        end
        if (s_z>0)
            l_w = sqrt((x-point_up(1))^2+(y-point_up(2))^2+(z-point_up(3))^2); %distance to lower boundary
        else 
            l_w = sqrt((x-point_down(1))^2+(y-point_down(2))^2+(z-point_down(3))^2); %distance to upper boundary
        end
        if point_1(3)==point_2(3)
            error('yatay')
        end


%         iceride_mi=sqrt(x*x+y*y+z*z)*10^6
%         l_w_micron=l_w*10^6
%         l_beta_micron = l_beta*10^6
%         if isnan(l_w)
%             error('dur')
%         end
        if l_w<l_beta %check if ray reaches to boundary or extinct?
            min_index=1; %reach boundary
            min_l=l_w;
        else
            min_index=2; %extinct (will absorbed or scatter? we will check below)
            min_l=l_beta;
            ever_extinct=1;
            donme_no=donme_no+1;
        end
% 
%         l_w_um=l_w*10^6
%         l_beta_um=l_beta*10^6

        x=x+min_l*s_x;% move
        y=y+min_l*s_y;% move
        z=z+min_l*s_z;% move

        if (min_index==1)
    %            disp('hit boundary');
            alive=0;
            if z>0
                t_no=1;%it left from bottom so transmitted
                if ever_extinct==0
                    t_direct_no=1;
                end
            else
                r_no=1;%it left from top so reflected
            end
        else
            random_no=rand();
            if random_no<scat_prob
    %               disp('scattering');
                if use_HG==1
                    [s_x,s_y,s_z]=scatter_hg(g,s_x,s_y,s_z);%find new trajectory by henyey greenstein phase function
                else
                    cos_theta= inv_cdf_cos(ceil(rand()*n_cdf_random)); %find the deviation from direction by scattering 
                    [s_x,s_y,s_z]=scatter_mc(cos_theta,s_x,s_y,s_z); %change the trajectory by cos_theta
                end
                l_beta=-log(rand())/mu_tot;%it extincted so don't keep the old l_beta create new for new event
            else
    %               disp('absorption');
                alive=0;%it is aborbed, game over for this ray. exit this loop
                a_no=1;
            end
        end
    end
    r_tot=r_tot+r_no;
    t_direct_tot=t_direct_tot+t_direct_no;
    t_tot=t_tot+t_no;
    a_tot=a_tot+a_no;
    donme_tot=donme_tot+donme_no;
    if ever_extinct==1
%         [theta,phi,z_] = cart2sph(x,y,z);
%         [theta,phi,z_] = cart2sph(s_x,s_y,s_z);
%         aci=round(phi*180/pi)+91;
%         aci=round(theta*180/pi)+181;
        aci=round(acosd(s_z));
        if aci==0
            aci=1;
        end
        acilar(aci)=acilar(aci)+1/sqrt(1-s_z*s_z);
    end
end
r_tot=r_tot/photon_number;
t_tot=t_tot/photon_number;
a_tot=a_tot/photon_number;
t_direct_tot=t_direct_tot/photon_number;
acilar=acilar/photon_number;
donme_sayisi=donme_tot/photon_number;
