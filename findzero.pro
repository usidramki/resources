


function rhs22, nu, lambda_B, rho_G_0, g_B_0, g_B

     oned0mnu = 1.d0 / (1.d0 - nu)
     twomnu = (2.d0 - nu) / (1.d0 - nu)

     kappa = (1.d0 - nu) * rho_G_0 * ( g_B_0 / (1.d0 - g_B_0) )^oned0mnu - nu * ibetam(g_B_0, oned0mnu, -oned0mnu) + lambda_B * ibetam(g_b_0, twomnu, -oned0mnu)

     rhs_eq = nu * ibetam(g_B, oned0mnu, -oned0mnu) - lambda_B * ibetam(g_B, twomnu, -oned0mnu) + kappa

     rho_G_eq = rhs_eq * (((1.d0 - g_B)/g_B)^oned0mnu) / (1.d0 - nu)

     func = rho_G_eq * (1.d0 - nu) * (1.d0 - g_B) * g_B

     return,func

;     return,rhs_eq

end




function solveforgB, nu, lambda_B, rho_G_0, g_B_0


npt_gB = 100

nzeros = 0
zeros = []

soln_x = dblarr(npt_gB)
soln_y = dblarr(npt_gB)

epsL=1d-4


for i_sc = 0, npt_gB-2 do begin

;        g_B_A = 1.d0-alog((1.d0*i_sc / npt_gB)+epsL)/alog(epsL)
;        g_B_B = 1.d0-alog(((i_sc+1.d0) / npt_gB)+epsL)/alog(epsL)

        g_B_A = 1.d0 * i_sc / npt_gB
        g_B_B = 1.d0 * (i_sc+1) / npt_gB



        rhs_vals0 = rhs22(nu, lambda_B, rho_G_0, g_B_0, g_B_A)
        rhs_vals1 = rhs22(nu, lambda_B, rho_G_0, g_B_0, g_B_B)
        
        soln_x[i_sc] = g_B_A
        soln_y[i_sc] = rhs_vals0

        if (rhs_vals0 * rhs_vals1) lt 0 then begin

;print,"QQQ=",rho_G_0,g_B_0,g_B_a,g_B_B,rhs_vals0, rhs_vals1

             nzeros++
             zeroval = g_B_A + (g_B_B - g_B_A) * rhs_vals0 / (rhs_vals1 - rhs_vals0)
             zeros = [zeros, zeroval]
        endif

endfor

plot,soln_x[0:-2],soln_y[0:-2],yrange=[-0.1,0.1]

as=''
read,as


return,[nzeros, zeros]

end










n_par_1 = 50
n_par_2 = 50


nu = 0.2d0
lambda_B = 0.3d0


n_zeros = dblarr(n_par_1, n_par_2)


for i_par_1 = 0, n_par_1-1 do begin
for i_par_2 = 0, n_par_2-1 do begin



rho_G_0 = 1.d0 * i_par_1 / n_par_1
g_B_0 = 1.d0 * i_par_2 / n_par_2


zeros = solveforgB(nu, lambda_B, rho_G_0, g_B_0)
print,zeros
n_zeros(i_par_1, i_par_2) = zeros[0]




endfor
print,i_par_1
endfor

tvframe,n_zeros,/ba,/as,xtitle='!7q!X!N!LG!X!N!E0!X!N',ytitle='g!LB!X!N!E0!X!N',charsize=3,xrange=[0,1],yrange=[0,1],title='!7m!X!N='+strcompress(string(nu,format='(f4.2)'),/remove_all)+', !7k!X!N!LB!X!N='+strcompress(string(lambda_B,format='(f4.2)'),/remove_all),btitle='# roots'


end

 
