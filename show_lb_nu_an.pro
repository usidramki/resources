pro openps,name
 SET_PLOT,'ps'
 DEVICE,FONT_SIZE=10
 nam_ps_file=name
 n_eps = 1
 sty='p'                ; p=portrait l=landscape
 phe=0.25 ; % of image height
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




device,decompose=0

restore,filename='zeros_down.sav'
restore,filename='lb_nu_scan_low_k.sav'


openps,'fig8a_s.eps'

!p.multi=[0,2,1]

loadct,0

tvframe2,BmR_scan_nu,xrange=[min(nu_grid),max(nu_grid)],yrange=[min(lambda_B_grid),max(lambda_B_grid)],xtitle='!7m!X!N',ytitle='!7k!X!N!LB!X!N',btitle='B-R',/bar,/sa,/as,charsize=cs

tek_color

contour,zeros_down,nu_grid,lambda_B_grid,/overplot,levels=[0.3],color=2,thick=5
contour,zeros_down,nu_grid,lambda_B_grid,/overplot,levels=[0],color=3,thick=5


loadct,0

tvframe,alog10(tte_scan_nu),xrange=[min(nu_grid),max(nu_grid)],yrange=[min(lambda_B_grid),max(lambda_B_grid)],xtitle='!7m!X!N',ytitle='!7k!X!N!LB!X!N',btitle='log(Time)',/bar,/sa,/as,charsize=cs

tek_color

contour,zeros_down,nu_grid,lambda_B_grid,/overplot,levels=[0.3],color=2,thick=5
contour,zeros_down,nu_grid,lambda_B_grid,/overplot,levels=[0],color=3,thick=5




closeps
spawn,'epspdf fig8a_s.eps fig8a_s.pdf'



end


