


function rhs22, nu, lambda_B, rho_G_0, g_B_0, g_B

     oned0mnu = 1.d0 / (1.d0 - nu)
     twomnu = (2.d0 - nu) / (1.d0 - nu)

     kappa = (1.d0 - nu) * rho_G_0 * ( g_B_0 / (1.d0 - g_B_0) )^oned0mnu - nu * ibetam(g_B_0, oned0mnu, -oned0mnu) + lambda_B * ibetam(g_b_0, twomnu, -oned0mnu)

     rhs_eq = nu * ibetam(g_B, oned0mnu, -oned0mnu) - lambda_B * ibetam(g_B, twomnu, -oned0mnu) + kappa

     rho_G_eq = rhs_eq * (((1.d0 - g_B)/g_B)^oned0mnu) / (1.d0 - nu)

;     func = rho_G_eq * (1.d0 - nu) * (1.d0 - g_B) * g_B

;     return,func

;     return,rhs_eq

return,rho_G_eq

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

        g_B_A = 1.d0 * i_sc / npt_gB; * 0.73
        g_B_B = 1.d0 * (i_sc+1) / npt_gB; * 0.73



        rhs_vals0 = rhs22(nu, lambda_B, rho_G_0, g_B_0, g_B_A)
        rhs_vals1 = rhs22(nu, lambda_B, rho_G_0, g_B_0, g_B_B)
        
        soln_x[i_sc] = g_B_A
        soln_y[i_sc] = rhs_vals0

        if (rhs_vals0 * rhs_vals1) lt 0 then begin

;print,"QQQ=",rho_G_0,g_B_0,g_B_a,g_B_B,rhs_vals0, rhs_vals1

             nzeros++
             zeroval = g_B_A + (g_B_B - g_B_A) * rhs_vals0 / (rhs_vals1 - rhs_vals0)

	     if (rhs_vals0 gt 0.d0) and (rhs_vals1 le 0.d0) then zeroval = -zeroval
	     
	     
	     
             zeros = [zeros, zeroval]
        endif

endfor

plot,soln_x[0:-2],soln_y[0:-2],xtitle='g_B', ytitle='rho_G',yrange=[-10,10]
;oplot,[soln_x[0],soln_x[-2]],[0,0],linestyle=2

;as = ''
;read,as

;stop

return,[nzeros, zeros]

end


npt_nu = 200
npt_lambda_B = 200

nu_grid = dindgen(npt_nu)/(npt_nu-1.d0)
lambda_B_grid = dindgen(npt_lambda_B)/(npt_lambda_B-1.d0)
kappa_B = 0.01d0
rho_G_0 = 0.1d0
g_B_0 = 0.3d0

zeros_down = dblarr(npt_nu, npt_lambda_B)
zeros_up = dblarr(npt_nu, npt_lambda_B)



for i_lambda_B = 0, npt_lambda_B-1 do begin
for i_nu = 0, npt_nu-1 do begin

	lambda_B = lambda_B_grid[i_lambda_B]
	nu = nu_grid[i_nu]

	zeros = solveforgB(nu, lambda_B, rho_G_0, g_B_0)
	print,zeros

        if (zeros[0] eq 0) then begin
	     zeros_down[i_nu, i_lambda_B] = -1
	     zeros_up[i_nu, i_lambda_B] = -1
	endif
    
	if (zeros[0] eq 1) then begin
	     ;if (zeros[1] lt 0) then zeros_down[i_nu, i_lambda_B] = abs(zeros[1]) else zeros_up[i_nu, i_lambda_B]=zeros[1]
	     zeros_down[i_nu,i_lambda_B] = abs(zeros[1]) 
	endif     	     
	     
;	if (zeros[0] ge 2) then begin
;	     if ((zeros[1] lt 0) and (zeros[2] gt 0)) then begin
;	             zeros_up[i_nu, i_lambda_B] = zeros[2]
;		     zeros_down[i_nu, i_lambda_B] = abs(zeros[1])
;             endif else begin
;	             zeros_up[i_nu, i_lambda_B] = zeros[1]
;		     zeros_down[i_nu, i_lambda_B] = abs(zeros[2])
;             endelse
;	endif		     

        if (zeros[0] eq 2) then zeros_down[i_nu,i_lambda_B] = abs(zeros[2])





endfor
endfor


tvframe,zeros_down,/ba,/as, xrange=[min(nu_grid),max(nu_grid)], yrange=[min(lambda_B_grid), max(lambda_B_grid)], xtitle='!7m!X!N', ytitle='!7k!X!N!LB!X!N'


save,filename='zeros_down.sav',lambda_b_grid,nu_grid,zeros_down



end

 
