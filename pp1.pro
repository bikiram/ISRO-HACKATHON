ff='newlxp2level2.evn'  
data=mrdfits(ff, 1)  
help, data, /struct  
time = data.TIME 
undefine,data
print,'time[0]',time[0]
;time=time[0:5000]
time=time-time[0]
print,'time[10]',time[10]
;end
;pro znscal,time,nulo,nuhi,dnu,iharm,ff,pw
nulo=0.10145d0
nuhi=0.10165d0
dnu=0.000001d0
iharm=8l
z=(nuhi-nulo)/dnu
print,'z',z
;end
pi2=2.0d0*!dpi
sum=n_elements(time)
print,'sum',sum
nnu=long((nuhi-nulo)/dnu)
print,'(nuhi-nulo)/dnu',(nuhi-nulo)/dnu
print,'nnu',nnu
;end
ff=dblarr(nnu)
pw=dblarr(nnu)
print,'total number of points = ',round(nnu)
print,'Beginning znscal with first of the trial freq.s'
j=0l
nu=0.0d0
for nu=nulo,nuhi-dnu,dnu do begin
	;print,'j1',j
        phi=pi2*nu*time
	print,n_elements(phi)
	;print,'phi',phi[0:10]
         zns=0.0
         for i=1l,iharm do begin
                a=total(cos(i*phi))

                b=total(sin(i*phi))
                zns=zns+(a*a+b*b)
	;print,'znsi',zns
         endfor
	print,'zns10',zns
	zns=zns*(2.0d0/sum)
	print,'zns',zns
        pw[j]=zns
        ff[j]=nu+0.5d0*dnu
        print,'j',j
	j=j+1l
endfor
print,'Completd znscal for each of the trial freq. ---> Proceeding to finding max. of them'
;print,'ff',ff
;print,'pw',pw
plot,ff,pw, xtitle='frequency(HZ)', ytitle='(Z^n)^2'
oplot,ff,pw,psym=1.0

!P.CHARSIZE=1 & !P.CHARTHICK=4
!P.THICK=3 & !X.THICK=6 & !Y.THICK=6
SET_PLOT,'PS'
;DEVICE,FILENAME='new_zns_obs2_har8.PS'
plot,ff,pw, xtitle='frequency(HZ)', ytitle='(Z^n)^2'
oplot,ff,pw,psym=2.0
al_legend,['OBS 2,k=8','nu_max=0.1015235'],psym=[1,1],box=2,/left,thick=2,spacing=2
DEVICE,/CLOSE
SET_PLOT,'X'
;y1 = histogram( pw,binsize=ffbin, locations = ff) 
;print,'y1',y1
;plot,ff,y1

end
