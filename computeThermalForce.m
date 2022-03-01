function [Fthermal_global,Fthermal_local,R_e,E_e,l_e,Ft] = computeThermalForce(Tmat,Temperature_increment,n_el,Td,n_d,x,n_dof,n_i,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%   - Temperature_increment
%   - n_el           Total number of elements
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fthermal_local:Force matrix caused by thermal expansion [4 x n_el]
%                   Treated as known forces
%                   Filas como fuerza en la dirección local [1x 1y 2x 2y]
%--------------------------------------------------------------------------
R_e = zeros(2*n_d, 2*n_d, n_el); %Sirve para pasar de base global a local
l_e = zeros(n_el,1);
E_e = zeros(n_el,1);
for i=1:n_el
    x_1_e= x(Tn(i,1),1);
    x_2_e= x(Tn(i,2),1);
    y_1_e= x(Tn(i,1),2);
    y_2_e= x(Tn(i,2),2);
    l_e(i)= sqrt((x_2_e-x_1_e)^2+(y_2_e-y_1_e)^2);
    s_e=(y_2_e-y_1_e)/l_e(i);
    c_e=(x_2_e-x_1_e)/l_e(i);
    mat_aux = [c_e s_e 0 0;
        -s_e c_e 0 0;
        0 0 c_e s_e;
        0 0 -s_e c_e
        ];
   R_e(:,:,i)= mat_aux;
   E_e(i) = Tmat(i, 1);
end

eps_thermal = zeros(n_el,1);
Fthermal_local = zeros(4,n_el);
for i = 1:n_el
    f= [ -1 0 1 0].'; %Vector de dirección de fuerzas local, únicamente en la dirección axial.
    eps_thermal(i,1)=Tmat(i,3)*Temperature_increment;
    Fthermal_local(:,i)=eps_thermal(i, 1)*Tmat(i,2)*Tmat(i,1)*f;
end
Fthermal_global = zeros(4,n_el);
for i = 1:n_el
    Fthermal_global(:,i) = (R_e(:,:,i))\Fthermal_local(:,i);
end
Ft = zeros(n_dof,1);
for i = 1:n_el
    for j = 1:(n_d*n_i)
    Ft((Td(i,j)),1)=Ft((Td(i,j)),1)+Fthermal_global(j,i);
    end
end

end