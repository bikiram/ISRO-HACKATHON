function modl,x,pars
;y=pars(0)+pars(1)*sin(2.0d0*!dpi*(x-pars(2)))
y=pars(0)+pars(1)*sin(2.0d0*!dpi*(x-pars(2)))+pars(3)*sin(4.0d0*!dpi*(x-pars(4)))+pars(5)*sin(6.0d0*!dpi*(x-pars(6)))+pars(7)*sin(8.0d0*!dpi*(x-pars(8)))+pars(9)*sin(10.0d0*!dpi*(x-pars(10)))+pars(11)*sin(12.0d0*!dpi*(x-pars(12)))+pars(13)*sin(14.0d0*!dpi*(x-pars(14)))+pars(15)*sin(16.0d0*!dpi*(x-pars(16)));+pars(17)*sin(18.0d0*!dpi*(x-pars(18)))+pars(19)*sin(20.0d0*!dpi*(x-pars(20)))+pars(21)*sin(22.0d0*!dpi*(x-pars(22)))+pars(23)*sin(24.0d0*!dpi*(x-pars(24))) +pars(25)*sin(26.0d0*!dpi*(x-pars(26)))+pars(27)*sin(28.0d0*!dpi*(x-pars(28)))+pars(29)*sin(30.0d0*!dpi*(x-pars(30))) ; harmonic
return, y
end

ff='newlxp2level2.evn'
undefine,data
data=mrdfits(ff, 1)
c1=52
c2=129
help, data, /struct 
time = data.TIME 
;print,tt
;end
ch=data.CHANNEL
time=time(where (ch ge c1 and ch lt c2))
help, data, /struct

pi2=2.0d0*!dpi
sum=n_elements(time)
print,time[0]
time=time-time[0]

nu2=0.1015235d0
phi1=nu2*time
undefine,data
print,'phi1',phi1[3000000:3000100]
phi2=phi1-long(phi1)

print,n_elements(phi2)
;print,'phi2',phi2[3000000:3000100]
;end
;phi3=alog10(phi2)
;print,'phi2',phi2[3000000:3000100]
;end
y1=dblarr(200l)
y1 = histogram( phi2, binsize = 0.001d0, locations = x, min=0.0d0,max=1.0d0)
z=n_elements(y1)
print,'z',z
phierr=sqrt(y1)
y1min=min(y1,iy1)
y2=shift(y1,-iy1)
phierr1=shift(phierr,-iy1)

;[code for creating phase folded light curve here]
result=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
;result = [700.0, -300.0,0.25,-300.0,0.25,300.0,-0.25,300.0,-0.25] ;initializes the fitting model parameters, can change this
result = mpfitfun('modl', x, y2, phierr1, result, bestnorm=chimin, perror=perrs, dof=dof)
result = mpfitfun('modl', x, y2, phierr1, result, bestnorm=chimin, perror=perrs, dof=dof)
result = mpfitfun('modl', x, y2, phierr1, result, bestnorm=chimin, perror=perrs, dof=dof)

print,'Parameters with error:'
print,result
print,perrs

print,'Reduced chi-square: ', chimin/(dof*1.0d0),'with dof of ', dof

amp=result[1]/result[0]
ampe= amp*sqrt( (perrs[1]/result[1])^2.0d0 + (perrs[0]/result[0])^2.0d0 )
print, 'Fractional amplitude = ', amp, '+/-', ampe

window,0,retain=2
plot,[x,1+x],[y2,y2],psym=10,/ys
errplot,[x,1+x],[y2-phierr1,y2-phierr1],[y2+phierr1,y2+phierr1]
oplot,[x,1+x],modl([x,1+x], result),thick=2,color=255
al_legend,['OBS 2, k=8'],box=0,/left,thick=2,spacing=2,color=mycolor('blue')
!P.CHARSIZE=1 & !P.CHARTHICK=4
!P.THICK=3 & !X.THICK=6 & !Y.THICK=6
SET_PLOT,'PS'
DEVICE,FILENAME='new_fitted_pp_obs2_har8_30_80kev.PS',/color
plot,[x,1+x],[y2,y2],psym=10,/ys, xtitle='phase', ytitle='count'
errplot,[x,1+x],[y2-phierr1,y2-phierr1],[y2+phierr1,y2+phierr1]
oplot,[x,1+x],modl([x,1+x], result),thick=2,color=mycolor('red')

al_legend,['OBS 2, k=8'],box=0,/left,thick=2,spacing=2,color=mycolor('blue')
DEVICE,/CLOSE
SET_PLOT,'X'
end
