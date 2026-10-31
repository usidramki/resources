pro openps,name
 SET_PLOT,'ps'
 DEVICE,FONT_SIZE=10
 nam_ps_file=name
 n_eps = 1
 sty='p'                ; p=portrait l=landscape
 phe=0.6 ; % of image height
 pwi=0.9   ;0.5   0.9   ; % if image width
 che=0.3      ; 0.3                ; % adjust height
 cwi=0.6
 psrat=[1.0,1.0]
 a_cpl='y'
 ps_dev, nam_ps_file, a_cpl, sty, phe,pwi,che,cwi, spn, ps_rat,eps=n_eps
end

function psi,x
alpha = 1d6
res = 0.5d0 * (tanh((x - 4.d0/alpha) * alpha) + 1.d0)
;if (x gt 1d-4) then res = x else res = 0.d0
return,res
end

function RHS, t, u, nu, lambda_B, lambda_R, kappa_B, kappa_R

;0 - rho_G_2d
;1 - g_B_2d

rho_G_2d = 0
g_B_2d = 1

omnu = 1.d0 - nu

rhs_vector = dblarr(n_elements(u))

rhs_vector[rho_G_2d] = u[rho_G_2d] * (nu - u[rho_G_2d] - lambda_B * u[g_B_2d])
rhs_vector[g_B_2d] = u[rho_G_2d] * omnu * (1.d0 - u[g_B_2d]) * u[g_B_2d]

return,rhs_vector
end












function solve_everything, nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0,t_1


         rho_G_2d = 0
         g_B_2d = 1

	 n_var = 2

	 initial_condition = dblarr(n_var) ; Delta_BG, Delta_GR, B, R, g, gamma


	 t_0 =  0.d0 ; beginning
         ;t_1 = 10.d0


	 initial_condition[rho_G_2d] = rho_G_0
         initial_condition[g_B_2d] = g_B_0
         
         

	 tol = 1.d-12 ;%%%

	 N_STAGE = 5



	 ;Fehlberg constants


	 if (N_STAGE eq 5) then begin

	 c = [0.d0, 1.d0/4.d0, 3.d0/8.d0, 12.d0/13.d0, 1.d0, 1.d0/2.d0]

	 a_rkf =  [[ 0.d0,              0.d0,            0.d0,             0.d0,             0.d0       ], $
        	   [ 1.D0/4.D0,         0.d0,            0.d0,             0.d0,             0.d0       ], $
        	   [ 3.D0/32.D0,        9.d0/32.d0,      0.d0,             0.d0,             0.d0       ], $
        	   [ 1932.D0/2197.D0,  -7200.d0/2197.d0, 7296.d0/2197.d0,  0.d0,             0.d0       ], $
        	   [ 439.D0/216.D0,    -8.d0,            3680.d0/513.d0,  -845.d0/4104.d0,   0.d0       ], $
        	   [-8.D0/27.D0,        2.d0,           -3544.d0/2565.d0,  1859.d0/4104.d0, -11.d0/40.d0]]




	 slope_weight =   [ 16.d0/135.d0, 0.d0, 6656.d0/12825.d0, 28561.d0/56430.d0, -9.d0/50.d0,  2.d0/55.d0]
	 err_weight =     [-1.d0 /360.d0, 0.d0, 128.d0/4275.d0,   2197.d0/75240.d0,  -1.d0/50.d0, -2.d0/55.d0]

	 endif


	 if (N_STAGE eq 6) then begin

	 c = [0.d0, 1.d0/5.d0, 3.d0/10.d0, 4.d0/5.d0, 8.d0/9.d0, 1.d0, 1.d0]

	 a_rkf = [[ 0.d0,              0.d0,            0.d0,             0.d0,             0.d0,             0.d0       ], $
        	  [ 1.d0/5.d0,         0.d0,            0.d0,             0.d0,             0.d0,             0.d0       ], $
        	  [ 3.d0/40.d0,        9.d0/40.d0,      0.d0,             0.d0,             0.d0,             0.d0       ], $
        	  [ 44.d0/45.d0,      -56.d0/15.d0,     32.d0/9.d0,       0.d0,             0.d0,             0.d0       ], $
        	  [ 19372.d0/6561.d0, -25360.d0/2187.d0,64448.d0/6561.d0,-212.d0/729.d0,    0.d0,             0.d0       ], $
        	  [ 9017.d0/3168.d0,  -355.d0/33.d0,    46732.d0/5247.d0, 49.d0/176.d0,    -5103.d0/18656.d0, 0.d0       ], $
		  [ 35.d0/384.d0,      0.d0,            500.d0/1113.d0,   125.d0/192.d0,   -2187.d0/6784.d0,  11.d0/84.d0]]


	 slope_weight = [35.d0/384.d0, 0.d0, 500.d0/1113.d0, 125.d0/192.d0, -2187.d0/6784.d0, 11.d0/84.d0, 0.d0]

	 err_weight = [71.d0/57600.d0, 0.d0, -71.d0/16695.d0, 71.d0/1920.d0, -17253.d0/339200.d0, 22.d0/525.d0, -1.d0/40.d0]

	 endif







	 ; Initialize
	 dt = (t_1 - t_0)
	 min_dt  = 1.d-10 ;%%% ;dt / 1d10
	 max_dt  = dt / 5d0;0.005;dt / 1d4
	 max_step = 30000L

	 tarr = t_0

	 dt = min_dt; dt / 100.d0
	 t = t_0
	 u = initial_condition
	 theta = u

	 du = dblarr(n_var)


	 need_to_make_step = 0L
	 n_singular = 0L
	 t_singular = [0]
	 i_step = 0L

	 slope = dblarr(n_elements(u), N_STAGE+1)

	 while ((t lt t_1) and (dt ge min_dt) and (i_step lt max_step)) do begin



		 if (t + dt gt t_1) then dt = t_1 - t

		 for i = 0,N_STAGE do begin

        	    for ii = 0,n_var-1 do begin
		      du[ii] = 0.d0
		      for jj = 0,N_STAGE-1 do begin
	        	  du[ii] = du[ii] + slope[ii,jj] * a_rkf[jj,i]
			 ; print,"a_rkf=",jj,i,a_rkf[jj,i]
        	      endfor
		    endfor

		    uuu = u + dt * du

		    ttt = t + c[i] * dt

		    slope[*, i] = RHS(ttt, uuu, nu, lambda_B, lambda_R, kappa_B, kappa_R)

        	 endfor



		 err_vec = dblarr(n_var)
		 for i = 0,n_var-1 do for jj = 0,N_STAGE do err_vec[i] = err_vec[i] + slope[i,jj] * err_weight[jj]


        	 err = sqrt( total ((dt*err_vec)^2.d0) )

		 if (finite(err) eq 0) then begin 
		    print,'Infinite'
		    break
		 endif

        	 err_bound = tol * max([sqrt(total(u^2.d0)), 1])




        	 if (err le err_bound) or (need_to_make_step eq 1) then begin
        	    t = t + dt
        	    u = u + dt * (slope # slope_weight)	   


        	    tarr   = [tarr,t]
        	    theta  = [[theta], [u]]

		    need_to_make_step = 0
        	    i_step++
		 endif

        	 if (err ne 0) then dt = min([max_dt, 0.8d0 * dt * (err_bound / err) ^ 0.2d0 ])

		 ;print, i_step, t, t_1, dt, min_dt, err, err_bound

        	 if (dt lt min_dt) then begin
		    print,'Singularity: forcing through. Precision lost.'
		    dt = min_dt
		    need_to_make_step = 1
		    n_singular++
		    t_singular = [t_singular, t]
		 endif

	 ; stop


         ; if (theta[B_8d_F,-1] lt 0.01) or (theta[R_8d_F,-1] lt 0.01) then break




	  endwhile

	  ;print, i_step,' N_singular = ',n_singular, ' times = ',t_singular


	  ;if ((i_step mod 100) eq 0) then if (t lt t_1) then print,"Either solver failed completely, or too many iterations."


          return, [reform(tarr,1,n_elements(tarr)), theta]

end







function rhs22, nu, lambda_B, rho_G_0, g_B_0, g_B

     oned0mnu = 1.d0 / (1.d0 - nu)
     twomnu = (2.d0 - nu) / (1.d0 - nu)

     kappa = (1.d0 - nu) * rho_G_0 * ( g_B_0 / (1.d0 - g_B_0) )^oned0mnu - nu * ibetam(g_B_0, oned0mnu, -oned0mnu) + lambda_B * ibetam(g_b_0, twomnu, -oned0mnu)

     rhs_eq = nu * ibetam(g_B, oned0mnu, -oned0mnu) - lambda_B * ibetam(g_B, twomnu, -oned0mnu) + kappa

     rho_G_eq = rhs_eq * (((1.d0 - g_B)/g_B)^oned0mnu) / (1.d0 - nu)


return,rho_G_eq

end





rho_G_2d = 0

g_B_2d = 1





rho_G_0 = 0.1d0
g_B_0 = 0.3d0


nu_grid = [0.1d0, 0.2d0, 0.5d0, 0.6d0]
lambda_B_grid = [0.05d0, 0.8d0, 0.6d0, 0.4d0]

rho_G_0 = 0.1d0
g_B_0 = 0.3d0

npt_lambda_B = n_elements(lambda_B_grid)
npt_nu = n_elements(nu_grid)

npt_gB = 20
grid_gB = dindgen(npt_gB) / (npt_gB - 1.d0)

npt_rhoG = 20
grid_rhoG = dindgen(npt_rhoG) / (npt_rhoG - 1.d0)

d_rho_G = dblarr(npt_gB, npt_rhoG)
d_g_B = dblarr(npt_gB, npt_rhoG)

val_gB = dblarr(npt_gB)


rho_G_max = 0.6

!p.multi=[0,3,4]

device,decompose=0
tek_color

openps,'fig9.eps'



cs1= 0.9
cs=1.2

for i_par= 0, npt_lambda_B-1 do begin

	lambda_B = lambda_B_grid[i_par]
	nu = nu_grid[i_par]

        omnu = 1.d0 - nu

        for i_gB = 0, npt_gB - 2 do val_gB[i_gB] = rhs22(nu, lambda_B, rho_G_0, g_B_0, grid_gB[i_gB])

        plot,grid_gB[0:-2], val_gB[0:-2], yrange=[-1,1],charsize=cs, xtitle='g!LB!X!N', ytitle='!7q!X!N!LG!X!N'

	
        xyouts,0.7,0.7,'!7m!X!N='+strcompress(string(nu,format='(f4.2)'),/remove_all),charsize=cs1
	xyouts,0.7,0.5,'!7k!X!N!LB!X!N='+strcompress(string(lambda_B,format='(f4.2)'),/remove_all),charsize=cs1



	qq1 = solve_everything(nu, lambda_B, 0.d0, 0.d0, 0.d0,  0.d0, 0.d0, rho_G_0, 0.d0, 0.d0, g_B_0, 0.d0, 100.d0)

	tarr1 = reform(qq1[0,*])
	theta1 = reform(qq1[1:*,*])

        plot,tarr1,theta1[rho_G_2d, *], xtitle='Time', ytitle='!7q!X!N!LG!X!N(t), g!LB!X!N(t)',/nodata,/xs,charsize=cs,yrange=[0,1]
	oplot,tarr1,theta1[rho_G_2d, *], color=3
        oplot,tarr1,theta1[g_B_2d, *], color=4


if (i_par eq 0) then begin
        xyouts,60,0.9,'!7q!X!N!LG!X!N',charsize=cs1
	oplot,[70,90],[0.9,0.9],color=3
	xyouts,60,0.78,'g!LB!X!N',charsize=cs1
	oplot,[70,90],[0.78,0.78],color=4
endif



        plot,[0,1],[0,rho_G_max],/nodata,charsize=cs, xtitle='g!LB!X!N', ytitle='!7q!X!N!LG!X!N'

        oplot, theta1[g_B_2d, *], theta1[rho_G_2d, *], thick=3

        oplot,grid_gB[0:-2], val_gB[0:-2], color=2, thick=3, linestyle=2

        for i_rhoG = 0, npt_rhoG-1 do begin
	   for i_gB = 0, npt_gB-1 do begin
	      rho_G_start = grid_rhoG[i_rhoG]
	      g_B_start = grid_gB[i_gB]
	      
	      if (rho_G_start lt rho_G_max) then begin
	      
	      d_rho_G[i_gB, i_rhoG] = rho_G_start * (nu - rho_G_start - lambda_B * g_B_start)
              d_g_B[i_gB, i_rhoG] = rho_G_start * omnu * (1.d0 -g_B_start) * g_B_start


kkk = 0.03
	      
	      vec_len = sqrt(d_rho_G[i_gB, i_rhoG]^2.d0 +  d_g_B[i_gB, i_rhoG]^2.d0)
	      d_rho_G[i_gB, i_rhoG] =  d_rho_G[i_gB, i_rhoG] / vec_len
	      d_g_B[i_gB, i_rhoG] =  d_g_B[i_gB, i_rhoG] / vec_len
	      
	      arrow, g_B_start, rho_G_start, g_B_start + d_g_B[i_gB, i_rhoG]*kkk, rho_G_start + d_rho_G[i_gB, i_rhoG]*kkk,/data,hsize=40	      


endif

	  	   
	   endfor
	endfor



endfor




closeps
spawn,'epspdf fig9.eps fig9_s.pdf'



end

 
