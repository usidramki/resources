pro openps,name
 SET_PLOT,'ps'
 DEVICE,FONT_SIZE=10
 nam_ps_file=name
 n_eps = 1
 sty='p'                ; p=portrait l=landscape
 phe=0.3 ; % of image height
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
;2 - rho_G_3d
;3 - g_B_3d
;4 - g_R_3d
;5 - B_8d_A
;6 - R_8d_A
;7 - rho_G_8d_A
;8 - rho_B_8d_A
;9 - rho_R_8d_A
;10 - g_G_8d_A
;11 - g_B_8d_A
;12 - g_R_8d_A
;13 - B_8d_F
;14 - R_8d_F
;15 - rho_G_8d_F
;16 - rho_B_8d_F
;17 - rho_R_8d_F
;18 - g_G_8d_F
;19 - g_B_8d_F
;20 - g_R_8d_F

rho_G_2d = 0
g_B_2d = 1
rho_G_3d = 2
g_B_3d = 3
g_R_3d = 4
B_8d_A = 5
R_8d_A = 6
rho_G_8d_A = 7
rho_B_8d_A = 8
rho_R_8d_A = 9
g_G_8d_A = 10
g_B_8d_A = 11
g_R_8d_A = 12
B_8d_F = 13
R_8d_F = 14
rho_G_8d_F = 15
rho_B_8d_F = 16
rho_R_8d_F = 17
g_G_8d_F = 18
g_B_8d_F = 19
g_R_8d_F = 20




omnu = 1.d0 - nu
kblb = kappa_B * lambda_B
krlr = kappa_R * lambda_R


rhs_vector = dblarr(n_elements(u))

rhs_vector[rho_G_2d] = u[rho_G_2d] * (nu - u[rho_G_2d] - lambda_B * u[g_B_2d])
rhs_vector[g_B_2d] = u[rho_G_2d] * omnu * (1.d0 - u[g_B_2d]) * u[g_B_2d]
rhs_vector[rho_G_3d] = (nu - u[rho_G_3d] - lambda_B * u[g_B_3d] - lambda_R * u[g_R_3d]) * u[rho_G_3d]
rhs_vector[g_B_3d] = u[rho_G_3d] * omnu * (1.d0 - u[g_B_3d] - u[g_R_3d]) * u[g_B_3d]
rhs_vector[g_R_3d] = (u[rho_G_3d] - nu) * (1.d0 - u[g_B_3d] - u[g_R_3d]) * u[g_R_3d]
rhs_vector[B_8d_A] = -u[rho_R_8d_A] * u[R_8d_A] * psi(u[B_8d_A])
rhs_vector[R_8d_A] = -u[rho_B_8d_A] * u[B_8d_A] * psi(u[R_8d_A])
rhs_vector[rho_G_8d_A] = (nu - u[rho_G_8d_A]) * u[rho_G_8d_A] - lambda_B * u[g_B_8d_A] * u[rho_G_8d_A] - lambda_R * u[g_R_8d_A] * u[rho_G_8d_A]
rhs_vector[rho_B_8d_A] = lambda_B * u[rho_G_8d_A] * u[g_B_8d_A] - u[B_8d_A] * u[rho_B_8d_A]
rhs_vector[rho_R_8d_A] = lambda_R * u[rho_G_8d_A] * u[g_R_8d_A] - u[R_8d_A] * u[rho_R_8d_A]
rhs_vector[g_B_8d_A] = (u[rho_G_8d_A] * u[B_8d_A] - nu * u[rho_B_8d_A]) * (1.d0 - u[g_B_8d_A] - u[g_R_8d_A]) * u[g_B_8d_A]
rhs_vector[g_R_8d_A] = (u[rho_G_8d_A] * u[R_8d_A] - nu) * (1.d0 - u[g_B_8d_A] - u[g_R_8d_A]) * u[g_R_8d_A]
rhs_vector[g_G_8d_A] = - rhs_vector[g_B_8d_A] - rhs_vector[g_R_8d_A]


rhs_vector[B_8d_F] = -u[rho_R_8d_F] * u[R_8d_F] * psi(u[B_8d_F])
rhs_vector[R_8d_F] = -u[rho_B_8d_F] * u[B_8d_F] * psi(u[R_8d_F])
rhs_vector[rho_G_8d_F] = (nu - u[rho_G_8d_F]) * u[rho_G_8d_F] - lambda_B * u[g_B_8d_F] * u[rho_G_8d_F] - lambda_R * u[g_R_8d_F] * u[rho_G_8d_F]
rhs_vector[rho_B_8d_F] = lambda_B * u[rho_G_8d_F] * u[g_B_8d_F] - u[B_8d_F] * u[rho_B_8d_F]
rhs_vector[rho_R_8d_F] = lambda_R * u[rho_G_8d_F] * u[g_R_8d_F] - u[R_8d_F] * u[rho_R_8d_F]

rhs_vector[g_B_8d_F] = u[g_B_8d_F] * (1.d0 - u[g_B_8d_F] - u[g_R_8d_F]) * (u[rho_G_8d_F] * u[B_8d_F] - nu * (1.d0 + kblb * (1.d0 - u[g_B_8d_F])) * u[rho_B_8d_F] / (1.d0 + kblb * (1.d0 - u[g_B_8d_F]) * u[rho_B_8d_F]))

rhs_vector[g_R_8d_F] = u[g_R_8d_F] * (1.d0 - u[g_B_8d_F] - u[g_R_8d_F]) * (u[rho_G_8d_F] * u[R_8d_F] - nu * (1.d0 + krlr * (1.d0 - u[g_R_8d_F])) * u[rho_R_8d_F] / (1.d0 + krlr * (1.d0 - u[g_R_8d_F]) * u[rho_R_8d_F]))

rhs_vector[g_G_8d_F] = - rhs_vector[g_B_8d_F] - rhs_vector[g_R_8d_F]
return,rhs_vector
end












function solve_everything, nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0, t_1


         rho_G_2d = 0
         g_B_2d = 1
         rho_G_3d = 2
         g_B_3d = 3
         g_R_3d = 4
         B_8d_A = 5
         R_8d_A = 6
         rho_G_8d_A = 7
         rho_B_8d_A = 8
         rho_R_8d_A = 9
         g_G_8d_A = 10
         g_B_8d_A = 11
         g_R_8d_A = 12
         B_8d_F = 13
         R_8d_F = 14
         rho_G_8d_F = 15
         rho_B_8d_F = 16
         rho_R_8d_F = 17
         g_G_8d_F = 18
         g_B_8d_F = 19
         g_R_8d_F = 20

	 n_var = 21

	 initial_condition = dblarr(n_var) ; Delta_BG, Delta_GR, B, R, g, gamma


	 t_0 =  0.d0 ; beginning



	 initial_condition[rho_G_2d] = rho_G_0
         initial_condition[g_B_2d] = g_B_0
         initial_condition[rho_G_3d] = rho_G_0
         initial_condition[g_B_3d] = g_B_0
         initial_condition[g_R_3d] = g_R_0
         initial_condition[B_8d_A] = B_0
         initial_condition[R_8d_A] = R_0
         initial_condition[rho_G_8d_A] = rho_G_0
         initial_condition[rho_B_8d_A] = rho_B_0
         initial_condition[rho_R_8d_A] = rho_R_0
         initial_condition[g_B_8d_A] = g_B_0
         initial_condition[g_R_8d_A] = g_R_0
         initial_condition[g_G_8d_A] = 1.d0 - g_B_0 - g_R_0
         initial_condition[B_8d_F] = B_0
         initial_condition[R_8d_F] = R_0
         initial_condition[rho_G_8d_F] = rho_G_0
         initial_condition[rho_B_8d_F] = rho_B_0
         initial_condition[rho_R_8d_F] = rho_R_0
         initial_condition[g_B_8d_F] = g_B_0
         initial_condition[g_R_8d_F] = g_R_0
         initial_condition[g_G_8d_F] = 1.d0 - g_B_0 - g_R_0         
         
         

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
	 min_dt  = 1.d-31 ;%%% ;dt / 1d10
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

		 print, i_step, t, t_1, dt, min_dt, err, err_bound

        	 if (dt lt min_dt) then begin
		    print,'Singularity: forcing through. Precision lost.'
		    dt = min_dt
		    need_to_make_step = 1
		    n_singular++
		    t_singular = [t_singular, t]
		 endif

	 ; stop

	  endwhile

	  print, i_step,' N_singular = ',n_singular, ' times = ',t_singular


	  if ((i_step mod 100) eq 0) then if (t lt t_1) then print,"Either solver failed completely, or too many iterations."


          return, [reform(tarr,1,n_elements(tarr)), theta]

end







device,decompose=0
tek_color





nu = 0.05d0
lambda_B = 0.05d0 
lambda_R = 1.d0 - lambda_B
kappa_B = 0.01d0
kappa_R = 1.d0 / kappa_B
B_0 = 1.d0
R_0 = 1.d0
rho_G_0 = 0.1d0
rho_B_0 = 0.0d0
rho_R_0 = 0.0d0
g_B_0 = 0.3d0
g_R_0 = 0.3d0






B_8d_A = 5
B_8d_F = 13

R_8d_A = 6
R_8d_F = 14

rho_G_2d = 0
rho_G_3d = 2
rho_G_8d_A = 7
rho_G_8d_F = 15

rho_B_8d_A = 8
rho_B_8d_F = 16

rho_R_8d_A = 9
rho_R_8d_F = 17

g_B_2d = 1
g_B_3d = 3
g_B_8d_A = 11
g_B_8d_F = 19

g_R_3d = 4
g_R_8d_A = 12
g_R_8d_F = 20


g_G_8d_A = 10
g_G_8d_F = 18









qq = solve_everything(nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0, 2000.d0)

tarr = reform(qq[0,*])
theta = reform(qq[1:*,*])







!p.multi=[0,1,3]






cs=3
th=1.5







plot,tarr, theta[B_8d_F, *], title='B R populations', xtitle='Time', ytitle='B(t), R(t)',/nodata,/xs,charsize=cs
;oplot,tarr, theta[B_8d_A, *], color=4, thick=th, linestyle=2
;oplot,tarr, theta[R_8d_A, *], color=2, thick=th, linestyle=2
oplot,tarr,theta[B_8d_F, *], color=4, thick=th
oplot,tarr,theta[R_8d_F, *], color=2, thick=th



plot,tarr,theta[rho_G_8d_F, *], title='Resources', xtitle='Time', ytitle='!7q!X!N!LGBR!X!N(t)',/nodata,/xs,charsize=cs
oplot,tarr,theta[rho_G_2d, *], color=3, thick=th,linestyle=1
;oplot,tarr,theta[rho_G_3d, *], color=3, thick=th,linestyle=3
;oplot,tarr,theta[rho_G_8d_A, *], color=3, thick=th,linestyle=2
oplot,tarr,theta[rho_G_8d_F, *], color=3, thick=th

;oplot,tarr,theta[rho_B_8d_A, *], color=4, thick=th, linestyle=2
oplot,tarr,theta[rho_B_8d_F, *], color=4, thick=th

;oplot,tarr,theta[rho_R_8d_A, *], color=2, thick=th, linestyle=2
oplot,tarr,theta[rho_R_8d_F, *], color=2, thick=th





plot,tarr,theta[g_B_8d_F, *], title='Green populations', xtitle='Time', ytitle='g!LGBR!X!N(t)',/nodata,/xs,charsize=cs
oplot,tarr,theta[g_B_2d, *], color=4, thick=th, linestyle=1
;oplot,tarr,theta[g_B_3d, *], color=4, thick=th, linestyle=3
;oplot,tarr,theta[g_B_8d_A, *], color=4, thick=th, linestyle=2
oplot,tarr,theta[g_B_8d_F, *], color=4, thick=th



;oplot,tarr,theta[g_R_3d, *], color=2, thick=th, linestyle=3
;oplot,tarr,theta[g_R_8d_A, *], color=2, thick=th, linestyle=2
oplot,tarr,theta[g_R_8d_F, *], color=2, thick=th


;oplot,tarr,theta[g_G_8d_A, *], color=3, thick=th, linestyle=2
oplot,tarr,theta[g_G_8d_F, *], color=3, thick=th


;window,1
;!p.multi=0

;rho_G_2d = 0
;rho_G_3d = 2
;rho_G_8d_A = 7
;rho_G_8d_F = 15


;plot,  theta[g_B_2d, *],  theta[rho_G_2d, *],   xrange=[0,1],yrange=[0,0.2],/nodata,xtitle='g_B', ytitle='rho_G'
;oplot, theta[g_B_2d, *],  theta[rho_G_2d, *],   thick=th, linestyle=1
;oplot, theta[g_B_3d, *],  theta[rho_G_3d, *],   thick=th, linestyle=3
;oplot, theta[g_B_8d_A, *],theta[rho_G_8d_A, *], thick=th, linestyle=4
;oplot, theta[g_B_8d_F, *],theta[rho_G_8d_F, *], thick=th






end
   
