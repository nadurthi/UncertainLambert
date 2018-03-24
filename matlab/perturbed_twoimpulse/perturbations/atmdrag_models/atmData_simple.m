function rho = atmData(h)
%%
% This function is based upon CIRA72_static   Static, Exponential Density Model
%   CIRA72_static calculates the local atmospheric density
%   using a static, exponentially decaying model. The static model is valid
%   for altitudes between 0 and 1000 km above the earth reference
%   ellipsoid. For altitutdes of 100 to 500 km the diurnal efect has been
%   averaged resulting in a mean daily density for this static model. 
%
%   References: Vallado - Fundamentals of Astrodynamics 3rd Edition, page
%                         562
%%

% Data    h (km)   \rho_0(kg/m^3)       H (km) % Table 8-4 from Vallado
Data =  [   0     1.225e0        8.44
           25     3.899e-2       6.49
           30     1.774e-2       6.75
           35     8.279e-3       7.07
           40     3.972e-3       7.47
           45     1.995e-3       7.83
           50     1.057e-3       7.95
           55     5.821e-4       7.73
           60     3.206e-4       7.29
           65     1.718e-4       6.81
           70     8.770e-5       6.33
           75     4.178e-5       6.00
           80     1.905e-5       5.70
           85     8.337e-6       5.41
           90     3.396e-6       5.38
           95     1.343e-6       5.74
          100     5.297e-7       6.15
          110     9.661e-8       8.06
          120     2.438e-8      11.6
          130     8.484e-9      16.1
          140     3.845e-9      20.6
          150     2.070e-9      24.6
          160     1.224e-9      26.3
          180     5.464e-10     33.2
          200     2.789e-10     38.5
          250     7.248e-11     46.9
          300     2.418e-11     52.5
          350     9.158e-12     56.4
          400     3.725e-12     59.4
          450     1.585e-12     62.2
          500     6.967e-13     65.8
          600     1.454e-13     79.0
          700     3.614e-14    109.0
          800     1.170e-14    164.0
          900     5.245e-15    225.0
         1000     3.019e-15    268.0];
     
rho = 0*h;
     
% Determine which altitude band object is in
for i =1:length(Data)
    dh = Data(i,1) - h;
    
    if dh == 0
        pointer = i;
        break
    elseif dh > 0
        pointer = i-1;
        break
    end
end

k=1;
% Determine Density
h0 = Data(pointer,1);
rho0 = Data(pointer,2);
H = Data(pointer,3);

if isequal(h0,dh) ==1
    rho(k) = rho0;
else
    rho(k) = rho0*exp((h0-h)/H);
end



end