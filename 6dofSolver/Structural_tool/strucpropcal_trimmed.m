clc
clear all
close all

format long g;

%-------------User Inputs------------------------------

base_dia = 0.214;

mass_warh = 91;

den_ext = 7800;

yield_stress = 850;

fos = 1.4;

int_pr = 14;

prop_load_fac = 0.8;

den_prop_opt = [1750 0];

no_fin = 4;

Nose_opt = 4;
    
Prop_opt = 1;   %(1 - Initial  2 - Burnout)

den_prop = den_prop_opt(Prop_opt);

%--------------- Thickness Calculation----------------------

tol = 1e-5;

err = 1;

max_stress = yield_stress/fos;

mean_dia = base_dia;

while (err > tol)

    t_noz_n = (int_pr * mean_dia) / (2 * max_stress); 

    int_dia = base_dia - 2 * t_noz_n;
        
    mean_dia = ( int_dia + base_dia )/2;
    
    t_noz = (int_pr * mean_dia) / (2 * max_stress);
    
    err = abs(t_noz - t_noz_n);
        
end

t_noz

if (base_dia <= 0.225)

    t_nose = 0.003; 
    
else
    
    t_nose = 0.004;
    
end

t_nose_ext = 0.005;

t_body = 0.003;

t_fin = 0.003; 

int_dia_nose = base_dia - 2 * t_nose;

int_dia_nose_ext = base_dia - 2 * t_nose_ext;

int_dia_body = base_dia - 2 * t_body;


%--------------- Warhead --------------------------

    %------------- Nose ---------------------------
    
    len_nose = 2.8 * base_dia; % factor
    
    switch (Nose_opt)
    
        case 1    %----------- Cone ---------------------- 

            ext_vol = (1/3) * pi * (base_dia/2)^2 * len_nose;

            int_vol = (1/3) * pi * (int_dia_nose/2)^2 * (len_nose - t_nose * 5); 
          
            cen_nose = (2/3) * len_nose ;
        
            iz_nose = (3/10) * (ext_vol - int_vol ) * den_ext * base_dia ^2;
    
            ix_nose = (3/80) * (ext_vol - int_vol ) * den_ext * (4*(base_dia/2)^2 + len_nose^2);
            
            
        case 2    %----------- Power Series - Parabola (n = 0.5) ---------------------- 

            ext_vol = (1/2) * pi * (base_dia/2)^2 * len_nose;

            int_vol = (1/2) * pi * (int_dia_nose/2 )^2 * (len_nose - t_nose * 3);
          
            cen_nose = (5/8) * len_nose;
            
            iz_nose = (3/10) * (ext_vol - int_vol ) * den_ext * base_dia ^2;  %%%%%%%%
    
            ix_nose = (3/80) * (ext_vol - int_vol ) * den_ext * (4*(base_dia/2)^2 + len_nose^2); %%%%%%%%%
            
            (ext_vol - int_vol ) * den_ext
        
        case 3    %----------- Elliptical ---------------------

            ext_vol = (2/3) * pi * (base_dia/2)^2 * len_nose;

            int_vol = (2/3) * pi * (int_dia_nose/2 )^2 * (len_nose - t_nose * 3);
          
            cen_nose =  len_nose - ((4/(3*pi)) * len_nose);

            iz_nose = (ext_vol - int_vol ) * den_ext * (base_dia/2)^2 / 5;
    
            ix_nose = (ext_vol - int_vol ) * den_ext * ((base_dia/2)^2 + len_nose^2)/ 10;
        
        case 4    %----------- Ogival ------------------------- 

            rad_n = ( (base_dia/2)^2 + len_nose^2 ) / base_dia;
        
            ext_vol = pi*(len_nose * rad_n ^2 - len_nose^3/3 - (rad_n - base_dia/2)* rad_n^2*asin(len_nose/rad_n));
        
            rad_n1 = ( (int_dia_nose/2)^2 + (len_nose-t_nose * 4)^2 ) / int_dia_nose;
        
            int_vol = pi*((len_nose-t_nose * 4) * rad_n1 ^2 - (len_nose-t_nose * 4)^3/3 - (rad_n1 - int_dia_nose/2)* rad_n1^2*asin((len_nose-t_nose * 4)/rad_n1));
        
            Mom_nose = 2 * pi * rad_n * ( 0.5 * len_nose^2 + (base_dia/2 - rad_n) * ( rad_n - sqrt(rad_n^2 - len_nose^2) ) );
        
            sa_nose = 2*pi*rad_n * ( ( (base_dia/2 - rad_n) * asin(len_nose/rad_n) ) + len_nose );
                
            cen_nose =  len_nose - (Mom_nose/sa_nose);
                
            iz_nose = (3/10) * (ext_vol - int_vol ) * den_ext * base_dia ^2;  %%%%%%%%
    
            ix_nose = (3/80) * (ext_vol - int_vol ) * den_ext * (4*(base_dia/2)^2 + len_nose^2); %%%%%%%%%
            
   
        case 5    %----------- Parabolic ---- ------------------
            
            ext_vol = (8/15) * pi * (base_dia/2)^2 * len_nose;

            int_vol = (8/15) * pi * (int_dia_nose/2 )^2 * (len_nose - t_nose * 4);
            
            cen_nose = (1/4) * len_nose;
        
            iz_nose = (3/10) * (ext_vol - int_vol ) * den_ext * base_dia ^2;  %%%%%%%%
    
            ix_nose = (3/80) * (ext_vol - int_vol ) * den_ext * (4*(base_dia/2)^2 + len_nose^2); %%%%%%%%%
            
        case 6   %-------------- Haack Series -----------------------
            
            C = (1/3) ;  % C = 0 for LD (von karman)  & C = 1/3 for LV
            
            ext_vol = simprl(@(x)haackfun(x,base_dia/2,len_nose,C),0,len_nose,10);
                      
            int_vol = simprl(@(x)haackfun(x,int_dia_nose/2,len_nose - t_nose * 4,C),0,len_nose-t_nose* 4,10);
            
            cen_nose = ((5/8) * len_nose + (1/4) * len_nose )/2 ; %%%%%%%%%%%
        
            iz_nose = (3/10) * (ext_vol - int_vol ) * den_ext * base_dia ^2;  %%%%%%%%
    
            ix_nose = (3/80) * (ext_vol - int_vol ) * den_ext * (4*(base_dia/2)^2 + len_nose^2); %%%%%%%%%
                    
    end
  
   
    mass_nose = (ext_vol - int_vol ) * den_ext
             
    %----------------Nose ext ----------------------

    len_nose_ext = 5.2 * base_dia; % factor

    ext_vol = pi * (base_dia/2)^2 * len_nose_ext;

    int_vol_1 = pi * (int_dia_nose_ext/2)^2 * len_nose_ext;

    mass_nose_ext = (ext_vol - int_vol_1 ) * den_ext;
          
    cen_nose_ext = (1/2) * len_nose_ext;
    
    mass_rat_nose = (int_vol / (int_vol_1 + int_vol) ) * (mass_warh - mass_nose * 2 - mass_nose_ext);
    
    mass_rat_nose_ext = (int_vol_1 / (int_vol_1 + int_vol) ) * (mass_warh - mass_nose - mass_nose_ext);

    
mass_warh_asm = mass_warh;
    
len_warh = len_nose + len_nose_ext;

iz_nose_ext = (mass_nose_ext + mass_rat_nose_ext) * (base_dia/2)^2 / 2;

ix_nose_ext = (1/12) * (mass_nose_ext + mass_rat_nose_ext) * (3 * (base_dia/2)^2 + len_nose_ext ^2);
    
%---------------- Body ---------------------------

len_body = 13.5 * base_dia; % factor

ext_vol = pi * (base_dia/2)^2 * len_body;

int_vol = pi * (int_dia_body/2)^2 * len_body;

mass_prop = den_prop * int_vol * prop_load_fac;

mass_ext = (ext_vol - int_vol )* den_ext  ;

mass_body_asm = mass_prop + mass_ext;

cen_body = (1/2) * len_body;

iz_body = (mass_body_asm) * (base_dia/2)^2 / 2;

ix_body = (1/12) * (mass_body_asm) * (3 * (base_dia/2)^2 + len_body ^2);

%------------- Nozzle --------------------------

len_noz = 1.8 * base_dia; % factor

ext_vol = pi * (base_dia/2)^2 * len_noz;

int_vol = pi * (int_dia/2)^2 * len_noz;

mass_noz = (ext_vol - int_vol )* den_ext * 2;

cen_noz = (1/2) * len_noz;

iz_noz = (mass_noz) * (base_dia/2)^2 / 2;

ix_noz = (1/12) * (mass_noz) * (3 * (base_dia/2)^2 + len_noz ^2);

%-------------- Fin ---------------------------

len_fin = 1.4 * base_dia; %factor

hei_fin = 0.7 * base_dia; %factor

vol_fin = len_fin * hei_fin * t_fin * no_fin * 1.1; % factor

mass_fin = vol_fin * den_ext;

cen_fin = (1/2) * len_fin;

iz_fin = (1/12) * (mass_fin) * (t_fin ^ 2 + hei_fin ^2) + mass_fin * (base_dia/2 + hei_fin/2)^2;

if ( no_fin == 0)
    
    no_fin = 1;
    
end

ix_fin = (1/12) * (mass_fin/no_fin) * 2 * (2 * len_fin^2 + t_fin^2 + hei_fin^2) + 2 * (mass_fin/no_fin) * (base_dia/2 + hei_fin/2)^2;

%------------- Assembly --------------------------

len_asm = len_warh + len_body + len_noz;

mass_asm = mass_warh_asm + mass_body_asm + mass_noz + mass_fin + 15 + mass_nose;

Total_Mass = mass_asm

cen_nose_ext_n = cen_nose_ext + len_nose;

cen_body_n = cen_body + len_warh;

cen_noz_n = cen_noz + len_warh + len_body;

cen_fin_n = cen_fin + len_warh + len_body;

cen_tot_asm = ( (cen_nose * (mass_nose * 2 + mass_rat_nose)) + (cen_nose_ext_n * (mass_nose_ext + mass_rat_nose_ext)) + (cen_body_n * mass_body_asm) + (cen_noz_n * mass_noz * 1.5) + (cen_fin_n * mass_fin * 2) ) / (mass_asm);

CG = cen_tot_asm

Ixx = ix_nose + (mass_nose * 2 + mass_rat_nose ) * (cen_tot_asm - cen_nose)^2 + ix_nose_ext + (mass_nose_ext + mass_rat_nose_ext) * (cen_tot_asm - cen_nose_ext_n)^2 + ix_body + (mass_body_asm) * (cen_tot_asm - cen_body_n)^2 + ix_noz + (mass_noz * 1.5) * (cen_tot_asm - cen_noz_n)^2 + ix_fin + (mass_fin * 2) * (cen_tot_asm - cen_fin_n)^2

Iyy = Ixx

Izz = iz_nose * 2 + iz_nose_ext + iz_body + iz_noz * 1.5 + iz_fin * 2


% ------------- Plots ----------------------------

% plot(cen_nose,base_dia,'*')
% hold on
% plot(cen_nose_ext_n , base_dia,'y*')
% hold on
% plot(cen_body_n ,base_dia,'r*')
% hold on
% plot(cen_noz_n ,base_dia,'k*')
% hold on
% plot(cen_fin_n ,base_dia,'m*')
% hold on
% plot(cen_tot_asm,base_dia,'g+')

%--------------------------------------------------







