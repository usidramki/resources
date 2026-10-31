function rhs22, nu, lambda_B, rho_G_0, g_B_0, g_B

     oned0mnu = 1.d0 / (1.d0 - nu)
     twomnu = (2.d0 - nu) / (1.d0 - nu)

     kappa = (1.d0 - nu) * rho_G_0 * ( g_B_0 / (1.d0 - g_B_0) )^oned0mnu - nu * ibetam(g_B_0, oned0mnu, -oned0mnu) + lambda_B * ibetam(g_b_0, twomnu, -oned0mnu)

     rhs_eq = nu * ibetam(g_B, oned0mnu, -oned0mnu) - lambda_B * ibetam(g_B, twomnu, -oned0mnu) + kappa

     rho_G_eq = rhs_eq / ((1.d0 - g_B)/g_B)^oned0mnu / (1.d0 - nu)

     func = rho_G_eq * (1.d0 - nu) * (1.d0 - g_B) * g_B

     return,func

end




function solveforgB, nu, lambda_B, rho_G_0, g_B_0




return,0

end










n_par_1 = 100
n_par_2 = 100


nu = 0.4d0
lambda_B = 0.3d0


rhs_vals = dblarr(n_par_1, n_par_2, 100)




for i_par_1 = 0, n_par_1-1 do begin
for i_par_2 = 0, n_par_2-1 do begin

for i_sc = 0, 99 do begin

g_B = 1.d0-alog((i_sc / 100.d0)+1d-3)/alog(1d-3)

rho_G_0 = 1.d0 * i_par_1 / n_par_1
g_B_0 = 1.d0 * i_par_2 / n_par_2


rhs_vals[i_par_1, i_par_2, i_sc] = rhs22(nu, lambda_B, rho_G_0, g_B_0, g_B)

;print,nu, lambda_B, rho_G_0, g_B_0, g_B, rhs_vals[i_par_1, i_par_2, i_sc]

endfor




endfor
print,i_par_1
endfor



end

 
