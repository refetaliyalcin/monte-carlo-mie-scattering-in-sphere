clc
clear all
close all


lambda=500*10^-9; %freespace wavelength of incident ray in meter 
R=100*10^-6;  %radius of medium in meter 
pigment_radius=150*10^-9; %radius of particles in meter 
f_v=0.01; %particle volume fraction. 0.01 corresponds to 1% 
use_HG=0; %if 0 use exact scattering phase function, if 1 uses henyey greenstein phase function approximation

n_pigment=1.5;  %real refractive index of pigment
k_pigment=0;  %imaginary refractive index of pigment
n_medium=1; %real refractive index of medium
k_medium=0; %imaginary refractive index of medium

photon_number=10^6; %number of rays that will be traced, higher the number more accurate the result
n_cdf_random=1000; %how many pieces will be between (0,1) in random number relation with inverse cumulative function, higher the number accurate the phase function. useless if HG is used
nang_gid=1000; %how many pieces will be between (0,pi) in cumulative distribution function, higher the number more accurate the result. useless if HG is used

lambda_nm=lambda*10^9; %for plot

%initialize variables
ref_lambda=zeros(length(lambda),1);
tra_lambda=zeros(length(lambda),1);
tra_direct_lambda=zeros(length(lambda),1);
abs_lambda=zeros(length(lambda),1);
acilar=zeros(180,length(lambda));

Qsca=zeros(length(lambda),1);
Qabs=zeros(length(lambda),1);
mu_tot=zeros(length(lambda),1);
scat_prob=zeros(length(lambda),1);
g=zeros(length(lambda),1);
if use_HG==0
    inv_cdf_cos=zeros(n_cdf_random,length(lambda));
    x2=flip(cos(transpose(linspace(0,pi,nang_gid))));
    zero_to_one=linspace(0,1,n_cdf_random)';
    zero_to_pi=linspace(0,pi,nang_gid)';
else
    inv_cdf_cos=zeros(1,length(lambda));
end

V=4*pi*pigment_radius^3/3; %volume of a single sphere
for i=1:length(lambda)
    x=2*pi*n_medium(i)*pigment_radius/lambda(i); %size parameter
    m=(n_pigment(i)+1i*k_pigment(i))/n_medium(i);
    fonksiyon=Mie(m,x); %calculate Lorenz Mie theory
    Qsca(i)=fonksiyon(2);%scattering efficiency
    Qabs(i)=fonksiyon(3);%absorption efficiency
    g(i)=fonksiyon(5);
    Csca=pi*pigment_radius^2*Qsca(i);%scattering crossection[m^2]
    Cabs=pi*pigment_radius^2*Qabs(i);%absorption crossection[m^2]
    alfa=f_v*Csca/V;%scattering coefficient[1/m]
    beta=f_v*Cabs/V+(1-f_v)*4*pi*k_medium(i)/lambda(i);%absorption coefficient[1/m] absorption of medium is implemented as a bulk property
    mu_tot(i)=alfa+beta;%extinction coefficient[1/m]
    scat_prob(i)=alfa/mu_tot(i);%scattering albedo also scattering probability
    if use_HG==0
        cdf=zeros(nang_gid,1);
        scat_ph_fn=Mie_phasefn(m, x, nang_gid)/(x*x*Qsca(i));%scattering phase function
        for i2=2:nang_gid
            cdf(i2)=trapz(x2(1:i2),scat_ph_fn(1:i2));%calculate cumulative distribution function
        end  
        inv_cdf_cos(:,i)=cos(interp1(cdf,zero_to_pi,zero_to_one,'linear','extrap'));%calculate inverse cumulative distribution function
    end
end
tic
%loop the monte carlo code for lambda and polar_angle
for i=1:length(lambda)
    [ref_lambda(i),tra_lambda(i),abs_lambda(i),tra_direct_lambda(i),acilar(:,i),donme(i)]=monte_carlo(photon_number,R,scat_prob(i),mu_tot(i),inv_cdf_cos(:,i),n_cdf_random,g(i),use_HG);
end
toc

scat_ph_fn=Mie_phasefn(m, x, 180)/(x*x*Qsca);%scattering phase function


set(0, 'DefaultLineLineWidth', 2); %set thickness of all the lines = 2
figure('Renderer', 'painters', 'Position', [500 300 428 420]) % starting point and height - width of the frame

set(gca, 'ColorOrder', [0 0 0;0 0 1; 1 0 0;0 0.5 0;0.75, 0.75, 0], 'NextPlot', 'replacechildren');% color of lines in the plot with the given order. remember it is periodic

hAx=gca;

degrees=linspace(0,pi,length(acilar));
% acilar(1:2)=acilar(3:4);
% acilar(end-1:end)=acilar(end-3:end-2);
acilar(isinf(acilar)|isnan(acilar)) = 0;
degerler=acilar/trapz(cos(degrees),-acilar);
plot(linspace(0,180,length(acilar)),degerler)
hold on
plot(linspace(0,180,length(scat_ph_fn)),scat_ph_fn)

y_max=1.05*max([max(scat_ph_fn),max(degerler)]);

hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
hLg=legend('Sphere','Single particle (Mie theory)','Location','east');
hLg.LineWidth=1.5;
hLg.EdgeColor = [0 0 0];
xlabel('Angle between entrance and exit [\circ]')
ylh=ylabel('Asymmetry factor, g');
ylh.VerticalAlignment	= 'bottom'; %if it is not alligned well, try 'top' and 'bottom' too
xlim([0 180])
ylim([0 y_max])
set(gca,'FontSize',13)
box on


g_mie=trapz(cos(degrees),-scat_ph_fn'.*cos(degrees));
g_mc=trapz(cos(degrees),-degerler'.*cos(degrees));
diffuse_tra=tra_lambda-tra_direct_lambda;
donme_compare=diffuse_tra+ref_lambda; %

averge_scattering_no=donme/donme_compare

txt1 = ['g_m_i_e = ' num2str(g_mie,3)];
txt2 = ['g_s_p_h = ' num2str(g_mc,3)];
txt3 = ['avg # of scat. event = ' num2str(averge_scattering_no,3)];
text(10,0.95*y_max,txt1,'FontSize',13)
text(10,0.85*y_max,txt2,'FontSize',13)
text(10,0.75*y_max,txt3,'FontSize',13)

xticks([0,45,90,135,180])
xticklabels({'0','45','90','135','180'})

saveas(gcf,'fig.png')
saveas(gcf,'fig.emf')
saveas(gcf,'fig.fig')




