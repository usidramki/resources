pro openps,name
 SET_PLOT,'ps'
 DEVICE,FONT_SIZE=10
 nam_ps_file=name
 n_eps = 1
 sty='p'                ; p=portrait l=landscape
 phe=0.53
 ; % of image height
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
rhs_vector[rho_G_8d_A] = (1.d0 - u[rho_G_8d_A]/nu) * u[rho_G_8d_A] - lambda_B * u[g_B_8d_A] * u[rho_G_8d_A] - lambda_R * u[g_R_8d_A] * u[rho_G_8d_A]
rhs_vector[rho_B_8d_A] = lambda_B * u[rho_G_8d_A] * u[g_B_8d_A] - u[B_8d_A] * u[rho_B_8d_A]
rhs_vector[rho_R_8d_A] = lambda_R * u[rho_G_8d_A] * u[g_R_8d_A] - u[R_8d_A] * u[rho_R_8d_A]
rhs_vector[g_B_8d_A] = (u[rho_G_8d_A] * u[B_8d_A] - nu * u[rho_B_8d_A]) * (1.d0 - u[g_B_8d_A] - u[g_R_8d_A]) * u[g_B_8d_A]
rhs_vector[g_R_8d_A] = (u[rho_G_8d_A] * u[R_8d_A] - nu) * (1.d0 - u[g_B_8d_A] - u[g_R_8d_A]) * u[g_R_8d_A]
rhs_vector[g_G_8d_A] = - rhs_vector[g_B_8d_A] - rhs_vector[g_R_8d_A]


rhs_vector[B_8d_F] = -u[rho_R_8d_F] * u[R_8d_F] * psi(u[B_8d_F])
rhs_vector[R_8d_F] = -u[rho_B_8d_F] * u[B_8d_F] * psi(u[R_8d_F])
rhs_vector[rho_G_8d_F] = (1.d0 - u[rho_G_8d_F]/nu) * u[rho_G_8d_F] - lambda_B * u[g_B_8d_F] * u[rho_G_8d_F] - lambda_R * u[g_R_8d_F] * u[rho_G_8d_F]
rhs_vector[rho_B_8d_F] = lambda_B * u[rho_G_8d_F] * u[g_B_8d_F] - u[B_8d_F] * u[rho_B_8d_F]
rhs_vector[rho_R_8d_F] = lambda_R * u[rho_G_8d_F] * u[g_R_8d_F] - u[R_8d_F] * u[rho_R_8d_F]

rhs_vector[g_B_8d_F] = u[g_B_8d_F] * (1.d0 - u[g_B_8d_F] - u[g_R_8d_F]) * (u[rho_G_8d_F] * u[B_8d_F] - nu * (1.d0 + kblb * (1.d0 - u[g_B_8d_F])) * u[rho_B_8d_F] / (1.d0 + kblb * (1.d0 - u[g_B_8d_F]) * u[rho_B_8d_F]))

rhs_vector[g_R_8d_F] = u[g_R_8d_F] * (1.d0 - u[g_B_8d_F] - u[g_R_8d_F]) * (u[rho_G_8d_F] * u[R_8d_F] - nu * (1.d0 + krlr * (1.d0 - u[g_R_8d_F])) * u[rho_R_8d_F] / (1.d0 + krlr * (1.d0 - u[g_R_8d_F]) * u[rho_R_8d_F]))

rhs_vector[g_G_8d_F] = - rhs_vector[g_B_8d_F] - rhs_vector[g_R_8d_F]

return,rhs_vector
end













function solve_everything, nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0,t_1

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
         ;t_1 = 10.d0


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


          if (theta[B_8d_F,-1] lt 0.01) or (theta[R_8d_F,-1] lt 0.01) or (theta[B_8d_F,-1] lt 0.01) or (theta[R_8d_F,-1] lt 0.01) then break




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


B_0 = 1.d0
R_0 = 1.d0
rho_B_0 = 0.1d0
rho_R_0 = 0.1d0
rho_G_0 = 0.1d0
g_R_0 = 0.3d0
g_B_0 = 0.3d0
g_G_0 = 1.d0 - g_R_0 - g_B_0


npt_lambda_B = n_elements(lambda_B_grid)
npt_nu = n_elements(nu_grid)

npt_gB = 200
grid_gB = dindgen(npt_gB) / (npt_gB - 1.d0)

npt_rhoG = 200
grid_rhoG = dindgen(npt_rhoG) / (npt_rhoG - 1.d0)

d_rho_G = dblarr(npt_gB, npt_rhoG)
d_g_B = dblarr(npt_gB, npt_rhoG)

val_gB = dblarr(npt_gB)


rho_G_max = 0.6

!p.multi=0

device,decompose=0
tek_color

openps,'fig11.eps'



cs1= 0.7
cs=1.2


npt_lambda_B = 200
npt_nu = 200
BmR_scan_nu_F = dblarr(npt_nu, npt_lambda_B)
BmR_scan_nu_A = dblarr(npt_nu, npt_lambda_B)

lambda_B_grid = dindgen(npt_lambda_B)/(npt_lambda_B-1.d0)
nu_grid = 10.d0*dindgen(npt_nu)/(npt_nu-1.d0)+1d-3

kappa_B = 0.01d0

calc=1
if (calc eq 1) then begin


     for i_lambda_B = 0, npt_lambda_B-1 do begin
     for i_nu = 0, npt_nu-1 do begin

     print,'solving for ',i_lambda_B, ' ', i_nu 

	     lambda_B = lambda_B_grid[i_lambda_B]
	     nu = nu_grid[i_nu]
	     lambda_R = 1.d0 - lambda_B
	     kappa_R = 1.d0 / kappa_B

	     qq1 = solve_everything(nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0, 5000.d0)
	     tarr1 = reform(qq1[0,*])
	     theta1 = reform(qq1[1:*,*])

             BmR_scan_nu_F(i_nu, i_lambda_B) = theta1[B_8d_F, -1] - theta1[R_8d_F, -1]

             BmR_scan_nu_A(i_nu, i_lambda_B) = theta1[B_8d_A, -1] - theta1[R_8d_A, -1]

	     compl = where((theta1[B_8d_F, *] ge 1d-2) * (theta1[R_8d_F, *] ge 1d-2) eq 0)

     endfor
     endfor

     save,filename='bmr_FA.sav', BmR_scan_nu_F, BmR_scan_nu_A, nu_grid, lambda_B_grid
endif else begin
      restore,filename='bmr_FA.sav'
endelse


dBmR = BmR_scan_nu_F - BmR_scan_nu_A
dBmR[0,*]=-0.001d0


tvframe2, dBmR,/ba,/sa,/as,btitle='[B-R]!UF!X!N-[B-R]!UA!X!N',xtitle='!7m!X!N',ytitle='!7k!X!N!LB!X!N',charsize=cs,xrange=[0,max(nu_grid)],yrange=[0,max(lambda_B_grid)]
tek_color
contour,dBmR,nu_grid,lambda_B_grid,levels=[-0.1],thick=2, color=0,/overplot,c_labels=[1]
contour,dBmR,nu_grid,lambda_B_grid,levels=[-0.05],thick=2, color=2,/overplot,c_labels=[1]
contour,dBmR,nu_grid,lambda_B_grid,levels=[0],thick=2, color=3,/overplot,c_labels=[1]
contour,dBmR,nu_grid,lambda_B_grid,levels=[0.05],thick=2, color=2,/overplot,c_labels=[1]
contour,dBmR,nu_grid,lambda_B_grid,levels=[0.1],thick=2, color=0,/overplot,c_labels=[1]






closeps
spawn,'epspdf fig11.eps fig11_s.pdf'



end

 
