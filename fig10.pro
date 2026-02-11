pro openps,name
 SET_PLOT,'ps'
 DEVICE,FONT_SIZE=10
 nam_ps_file=name
 n_eps = 1
 sty='p'                ; p=portrait l=landscape
 phe=0.2 ; % of image height
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

rhs_vector[B_8d_A] = -u[rho_R_8d_A] * u[R_8d_A] ;* psi(u[B_8d_A])
rhs_vector[R_8d_A] = -u[rho_B_8d_A] * u[B_8d_A] ;* psi(u[R_8d_A])
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


         ; if (theta[B_8d_F,-1] lt 0.01) or (theta[R_8d_F,-1] lt 0.01) then break




	  endwhile

	  ;print, i_step,' N_singular = ',n_singular, ' times = ',t_singular


	  ;if ((i_step mod 100) eq 0) then if (t lt t_1) then print,"Either solver failed completely, or too many iterations."


          return, [reform(tarr,1,n_elements(tarr)), theta]

end







;#############################################################################
;
; Copyright (C) 2008, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment and a link to the website.
; 
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;
; hf = Hypergeometric2F1(a, b, c, x)
;
; From the book "Computation of Special Functions"
;    by Shanjie Zhang and Jianming Jin
;    Copyright 1996 by John Wiley & Sons, Inc.
; The authors state:
;   "We give permission to the reader who purchases this book 
;    to incorporate any of these programs into his or her 
;    programs provided that the copyright is acknowledged."
;
;       ====================================================
;       Purpose: Compute hypergeometric function F(a,b,c,x)
;       Input :  a --- Parameter
;                b --- Parameter
;                c --- Parameter, c <> 0,-1,-2,...
;                x --- Argument   ( x < 1 )
;       Output:  HF --- F(a,b,c,x)
;       Routines called:
;            (1) GAMMA for computing gamma function
;            (2) PSI for computing psi function
;       ====================================================
;
;------------------------------------------------------------------------------
function hygfx_psi, x
compile_opt idl2, hidden
;
;       ======================================
;       Purpose: Compute PolyGamma Psi function
;       Input :  x  --- Argument of psi(x)
;       Output:  PS --- psi(x)
;       ======================================

el = 0.57721566490153286061d ; N[EulerGamma,20] in Mathematica

xa = ABS(x)
s = 0.0D0
IF (x eq fix(x) and x LE 0.0) THEN begin
  ps = 1.0D+300
  RETURN, ps
endif ELSE IF (xa eq fix(xa)) THEN begin
  n = xa
  for  k = 1, n - 1 do s = s + 1.0D0 / k
  ps = -el + s
endif ELSE IF (xa+0.5 eq fix(xa+0.5)) THEN begin
  n = xa - 0.5d 
  for  k = 1, n do s = s + 1.0d / (2*k-1)
  ps = -el + 2.0D0 * s - 1.386294361119891D0
endif else begin
  IF (xa LT 10.0) THEN begin
    n = 10 - fix(xa)
    for  k = 0, n - 1 do s = s + 1.0D0 / (xa + k)
    xa = xa + n
  endif
  x2 = 1.0D0 / (xa*xa)
  a1 = -.83333333333333333D-01
  a2 =  .83333333333333333D-02
  a3 = -.39682539682539683D-02
  a4 =  .41666666666666667D-02
  a5 = -.75757575757575758D-02
  a6 =  .21092796092796093D-01
  a7 = -.83333333333333333D-01
  a8 =  .4432598039215686D0
  ps = aLOG(xa) - .5D0 / xa + x2 * (((((((a8*x2 + a7)*x2 + a6)*x2 + a5)*  $
       x2 + a4)*x2 + a3)*x2 + a2)*x2 + a1)
  ps = ps - s
endelse
IF (x LT 0.0) then ps = ps - pi * COS(pi*x) / SIN(pi*x) - 1.0D0 / x
RETURN, ps

END
;------------------------------------------------------------------------------
function hygfx_Hypergeometric2F1, a, b, c, x
compile_opt idl2
on_error, 2

el = 0.57721566490153286061d ; N[EulerGamma,20] in Mathematica

l0 = c eq fix(c) and c LT 0.0
l1 = 1.0D0 - x LT 1.0D-15 and c - a - b LE 0.0
l2 = a eq fix(a) and a LT 0.0
l3 = b eq fix(b) and b LT 0.0
l4 = c - a eq fix(c-a) and c - a LE 0.0
l5 = c - b eq fix(c-b) and c - b LE 0.0
IF (l0 or l1) THEN message, 'The hypergeometric series is divergent'
eps = 1.0D-15
IF (x GT 0.95) then eps = 1.0D-8
IF (x eq 0.0 or a eq 0.0 OR b eq 0.0) THEN begin
  hf = 1.0D0
  return, hf
endif ELSE IF (1.0D0-x eq eps and c-a-b GT 0.0) THEN begin
  gc = gamma(c)
  gcab = gamma(c-a-b)
  gca = gamma(c-a)
  gcb = gamma(c-b)
  hf = gc * gcab / (gca*gcb)
  return, hf
endif ELSE IF (1.0D0+x LE eps and ABS(c-a+b-1.0) LE eps) THEN begin
  g0 = SQRT(!dpi) * 2.0D0^(-a)
  g1 = gamma(c)
  g2 = gamma(1.0D0+a/2.0-b)
  g3 = gamma(0.5D0+0.5*a)
  hf = g0 * g1 / (g2*g3)
  return, hf
endif ELSE IF (l2 or l3) THEN begin
  IF (l2) then nm = fix(ABS(a))
  IF (l3) then nm = fix(ABS(b))
  hf = 1.0D0
  r = 1.0D0
  for  k = 1, nm do begin
    r = r * (a+k-1.0D0) * (b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
  endfor
  return, hf
endif ELSE IF (l4 or l5) THEN begin
  IF (l4) then nm = fix(ABS(c-a))
  IF (l5) then nm = fix(ABS(c-b))
  hf = 1.0D0
  r = 1.0D0
  for  k = 1, nm do begin
    r = r * (c-a+k-1.0D0) * (c-b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
  endfor
  hf = (1.0D0-x)^(c-a-b) * hf
  return, hf
endif
aa = a
bb = b
x1 = x
IF (x LT 0.0D0) THEN begin
  x = x / (x-1.0D0)
  IF (c GT a and b LT a and b GT 0.0) THEN begin
    a = bb
    b = aa
  endif
  b = c - b
endif
IF (x GE 0.75D0) THEN begin
  gm = 0.0D0
  IF (ABS(c-a-b-fix(c-a-b)) LT 1.0D-15) THEN begin
    m = fix(c-a-b)
    ga = gamma(a)
    gb = gamma(b)
    gc = gamma(c)
    gam = gamma(a+m)
    gbm = gamma(b+m)
    pa = hygfx_psi(a)
    pb = hygfx_psi(b)
    IF (m ne 0) then gm = 1.0D0
    for  j = 1, ABS(m) - 1 do gm = gm * j
    rm = 1.0D0
    for  j = 1, ABS(m) do rm = rm * j
    f0 = 1.0D0
    r0 = 1.0D0
    r1 = 1.0D0
    sp0 = 0.0d0
    sp = 0.0D0
    IF (m ge 0) THEN begin
      c0 = gm * gc / (gam*gbm)
      c1 = -gc * (x-1.0D0)^m / (ga*gb*rm)
      for  k = 1, m - 1 do begin
        r0 = r0 * (a+k-1.0D0) * (b+k-1.0) / (k*(k-m)) * (1.0-x)
        f0 = f0 + r0
      endfor
      for  k = 1, m do $
        sp0 = sp0 + 1.0D0 / (a+k-1.0) + 1.0 / (b+k-1.0) - 1.0 / k
      f1 = pa + pb + sp0 + 2.0D0 * el + aLOG(1.0D0-x)
      hw = f1
      for  k = 1, 250 do begin
        sp = sp + (1.0D0-a) / (k*(a+k-1.0)) + (1.0-b) / (k*(b+k-1.0))
        sm = 0.0D0
        for  j = 1, m do $
          sm = sm + (1.0D0-a) / ((j+k)*(a+j+k-1.0)) + 1.0 / (b+j+k-1.0)
        rp = pa + pb + 2.0D0 * el + sp + sm + aLOG(1.0D0-x)
        r1 = r1 * (a+m+k-1.0D0) * (b+m+k-1.0) / (k*(m+k)) * (1.0-x)
        f1 = f1 + r1 * rp
        IF (ABS(f1-hw) LT ABS(f1)*eps) then break 
        hw = f1
      endfor
      hf = f0 * c0 + f1 * c1
    endif ELSE IF (m LT 0) THEN begin
      m = -m
      c0 = gm * gc / (ga*gb*(1.0D0-x)^m)
      c1 = -(-1)^m * gc / (gam*gbm*rm)
      for  k = 1, m - 1 do begin
        r0 = r0 * (a-m+k-1.0D0) * (b-m+k-1.0) / (k*(k-m)) * (1.0-x)
        f0 = f0 + r0
      endfor
      for  k = 1, m do sp0 = sp0 + 1.0D0 / k
      f1 = pa + pb - sp0 + 2.0D0 * el + aLOG(1.0D0-x)
      hw = f1
      for  k = 1, 250 do begin
        sp = sp + (1.0D0-a) / (k*(a+k-1.0)) + (1.0-b) / (k*(b+k-1.0))
        sm = 0.0D0
        for  j = 1, m do sm = sm + 1.0D0 / (j+k)
        rp = pa + pb + 2.0D0 * el + sp - sm + aLOG(1.0D0-x)
        r1 = r1 * (a+k-1.0D0) * (b+k-1.0) / (k*(m+k)) * (1.0-x)
        f1 = f1 + r1 * rp
        IF (ABS(f1-hw) LT ABS(f1)*eps) then break 
        hw = f1
      endfor
      hf = f0 * c0 + f1 * c1
    endif
  endif ELSE begin
    ga = gamma(a)
    gb = gamma(b)
    gc = gamma(c)
    gca = gamma(c-a)
    gcb = gamma(c-b)
    gcab = gamma(c-a-b)
    gabc = gamma(a+b-c)
    c0 = gc * gcab / (gca*gcb)
    c1 = gc * gabc / (ga*gb) * (1.0D0-x)^(c-a-b)
    hf = 0.0D0
    hw = hf
    r0 = c0
    r1 = c1
    for  k = 1, 250 do begin
      r0 = r0 * (a+k-1.0D0) * (b+k-1.0) / (k*(a+b-c+k)) * (1.0-x)
      r1 = r1 * (c-a+k-1.0D0) * (c-b+k-1.0) / (k*(c-a-b+k)) * (1.0-x)
      hf = hf + r0 + r1
      IF (ABS(hf-hw) LT ABS(hf)*eps) then break
      hw = hf
    endfor
    hf = hf + c0 + c1
  endelse
endif ELSE begin
  a0 = 1.0D0
  IF (c GT a and c LT 2.0D0*a and c GT b and c LT 2.0D0*b) THEN begin
    a0 = (1.0D0-x)^(c-a-b)
    a = c - a
    b = c - b
  endif
  hf = 1.0D0
  hw = hf
  r = 1.0D0
  for  k = 1, 500 do begin
    r = r * (a+k-1.0D0) * (b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
    IF (ABS(hf-hw) LE ABS(hf)*eps) then break 
    hw = hf
  endfor
  hf = a0 * hf
endelse
IF (x1 LT 0.0D0) THEN begin
  x = x1
  c0 = 1.0D0 / (1.0D0-x)^aa
  hf = c0 * hf
endif
a = aa
b = bb
;IF (k GT 120) then print, ' Warning; You should check the accuracy'
return, hf

END
;------------------------------------------------------------------------------
function Hypergeometric2F1, a, b, c, x
compile_opt idl2
on_error, 2
;
; This is just a wrapper for the actual function, 
; as the computation routine is not (yet) vectorized.

n = n_elements(x)
hyp = dblarr(n,/NOZERO)
for j=0,n-1 do hyp[j] = hygfx_Hypergeometric2F1(a,b,c,x[j])

return, hyp
END

 
 
;----------------------------------------------------------------------
function ibetam, x, a, b
;
; Incomplete beta function defined in the same 
; way as the Mathematica function Beta[x,a,b].
; This routine works for x<1 and was tested against Mathematica.
;
; V1.0: Michele Cappellari, Oxford, 01/APR/2008
; V2.0: Use Hypergeometric function for negative parameters.
;    From equation (6.6.8) of Abramoviz & Stegun
;    http://www.nrbook.com/abramowitz_and_stegun/page_263.htm 
;    MC, Oxford, 04/APR/2008
 
if a gt 0 && b gt 0 then $
    ib = ibeta(a,b,x)*beta(a,b) $
else $ 
    ib = x^a/a * hypergeometric2F1(a, 1d - b, a + 1d, x)
return, ib 
end
;----------------------------------------------------------------------



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

npt_gB = 40
grid_gB = dindgen(npt_gB) / (npt_gB - 1.d0)

npt_rhoG = 80
grid_rhoG = dindgen(npt_rhoG) / (npt_rhoG - 1.d0)

d_rho_G = dblarr(npt_gB, npt_rhoG)
d_g_B = dblarr(npt_gB, npt_rhoG)

val_gB = dblarr(npt_gB)


rho_G_max = 0.6

!p.multi=[0,3,1]

device,decompose=0
tek_color

;openps,'fig10.eps'
openps,'fig10_test.eps'


cs1= 0.7
cs=1.2

	lambda_B = 0.9d0
	lambda_R = 1.d0 - lambda_B
	nu = 0.2d0
	
	kappa_B = 0.001d0
	kappa_R = 1.d0 / kappa_B

        omnu = 1.d0 - nu



	qq1 = solve_everything(nu, lambda_B, lambda_R, kappa_B, kappa_R,  B_0, R_0, rho_G_0, rho_B_0, rho_R_0, g_B_0, g_R_0, 50.d0)

	tarr1 = reform(qq1[0,*])
	theta1 = reform(qq1[1:*,*])

        plot,tarr1,theta1[rho_G_2d, *], xtitle='Time', ytitle='!7q!X!N!LG!X!N(t), g!LB!X!N(t)',/nodata,/xs,charsize=cs,yrange=[0,1],title='2D reduced system'
	oplot,tarr1,theta1[rho_G_2d, *], color=3
        oplot,tarr1,theta1[g_B_2d, *], color=4



        xyouts,35,0.6,'!7q!X!N!LG!X!N',charsize=cs1
	oplot,[40,45],[0.6,0.6],color=3
	xyouts,35,0.55,'g!LB!X!N',charsize=cs1
	oplot,[40,45],[0.55,0.55],color=4



plot,tarr1, theta1[B_8d_F, *], title='Full system', xtitle='Time', ytitle='B(t), R(t)',/nodata,/xs,charsize=cs,yrange=[0,1]
oplot,tarr1,theta1[B_8d_F, *], color=4, thick=th, linestyle = 1
oplot,tarr1,theta1[R_8d_F, *], color=2, thick=th, linestyle = 1
oplot,tarr1,theta1[rho_G_8d_F, *], color=3, thick=th
oplot,tarr1,theta1[g_B_8d_F, *], color=4, thick=th

oplot,tarr1,theta1[B_8d_A, *], color=4, thick=th, linestyle = 2
oplot,tarr1,theta1[R_8d_A, *], color=2, thick=th, linestyle = 2
oplot,tarr1,theta1[rho_G_8d_A, *], color=3, thick=th, linestyle = 2
oplot,tarr1,theta1[g_B_8d_A, *], color=4, thick=th, linestyle = 3
oplot,tarr1,1.d0-theta1[g_B_8d_A, *]-theta1[g_R_8d_A, *], color=3, thick=th, linestyle = 4
;oplot,tarr1,theta1[g_R_8d_A,*],thick=th
;oplot,tarr1,theta1[g_R_8d_F,*],thick=th,linestyle=2



qq = 0.01+0.3
qq1 = 0.3

        xyouts,15,0.6-qq,'B!UF!X!N',charsize=cs1
	oplot,[20,25],[0.6-qq1,0.6-qq1],color=4, thick=th, linestyle = 1
        xyouts,15,0.55-qq,'B!UA!X!N',charsize=cs1
	oplot,[20,25],[0.55-qq1,0.55-qq1],color=4, thick=th, linestyle = 2
        xyouts,15,0.5-qq,'R!UF!X!N',charsize=cs1
	oplot,[20,25],[0.5-qq1,0.5-qq1], color=2, thick=th, linestyle = 1
        xyouts,15,0.45-qq,'R!UA!X!N',charsize=cs1
	oplot,[20,25],[0.45-qq1,0.45-qq1], color=2, thick=th, linestyle = 2
		
	xyouts,15,0.4-qq,'g!LB!UF!X!N',charsize=cs1
	oplot,[20,25],[0.4-qq1,0.4-qq1], color=4, thick=th


qq = 0.01

	xyouts,35,0.3-qq,'g!LB!UA!X!N',charsize=cs1
	oplot,[40,45],[0.3,0.3], color=4, thick=th, linestyle = 3

        xyouts,35,0.2-qq,'!7q!X!N!LG!UF!X!N',charsize=cs1
	oplot,[40,45],[0.2,0.2], color=3, thick=th
        xyouts,35,0.1-qq,'!7q!X!N!LG!UA!X!N',charsize=cs1
	oplot,[40,45],[0.1,0.1], color=3, thick=th, linestyle = 2








;        plot,[0,1],[0,rho_G_max],/nodata,charsize=cs, xtitle='g!LB!X!N', ytitle='!7q!X!N!LG!X!N'

        plot,[0,0.6],[0,0.2],/nodata,charsize=cs, xtitle='g!LB!X!N', ytitle='!7q!X!N!LG!X!N'


        oplot, theta1[g_B_2d, *], theta1[rho_G_2d, *], thick=3

        ;oplot,grid_gB[0:-2], val_gB[0:-2], color=2, thick=3, linestyle=2

        for i_rhoG = 0, npt_rhoG-1 do begin
	   for i_gB = 0, npt_gB-1 do begin
	      rho_G_start = grid_rhoG[i_rhoG]
	      g_B_start = grid_gB[i_gB]
	      
	      if (rho_G_start lt rho_G_max) then begin

		  d_rho_G[i_gB, i_rhoG] = rho_G_start * (nu - rho_G_start - lambda_B * g_B_start)
        	  d_g_B[i_gB, i_rhoG] = rho_G_start * omnu * (1.d0 -g_B_start) * g_B_start

        	  kkk = 0.01

		  vec_len = sqrt(d_rho_G[i_gB, i_rhoG]^2.d0 +  d_g_B[i_gB, i_rhoG]^2.d0)
		  d_rho_G[i_gB, i_rhoG] =  d_rho_G[i_gB, i_rhoG] / vec_len
		  d_g_B[i_gB, i_rhoG] =  d_g_B[i_gB, i_rhoG] / vec_len

		  if (g_B_start lt 0.6) and (rho_G_start lt 0.2) then arrow, g_B_start, rho_G_start, g_B_start + d_g_B[i_gB, i_rhoG]*kkk, rho_G_start + d_rho_G[i_gB, i_rhoG]*kkk,/data,hsize=60
	      
	      endif
	  	   
	   endfor
	endfor











closeps
spawn,'epspdf fig10_test.eps fig10_s_test.pdf'



end

 
