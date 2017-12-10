
pro spiralpot


;-----set up window size, window title-----

window, 0, xsize=800, ysize=800, title='polar-contour'


;-----set device & colors-----

device, decompose=0
device, true_color=24
;!p.color=0
;!p.background=16777215
;!p.background=255


;---number of grid----

nx=240
ny=64


;-----r, phi, rho array definition-----

r=fltarr(nx)
phi=fltarr(ny)
spiralpot=fltarr(ny, nx)


;-----reading coordinates-----

openr, lun, './idl0/r.txt', /get_lun
readf, lun, r

r=r*10


openr, lun, './idl0/phi.txt', /get_lun
readf, lun, phi



;-----choosing color table-----

;loadct, 0    
;loadct, 1
;loadct, 2    
;loadct, 3
;loadct, 4
;loadct, 5
;loadct, 6    
;loadct, 7
;loadct, 8
;loadct, 10
;loadct, 11
;loadct, 12
;loadct, 13
;loadct, 14
;loadct, 15
;loadct, 16
loadct, 17
;loadct, 18 ;only 41
;loadct, 22
;loadct, 23
;loadct, 24
;loadct, 25
;loadct, 33
;loadct, 34
;loadct, 38
;loadct, 39
;loadct, 3, file='/Applications/exelis/idl/resource/colors/colors2.tbl'


;-----set value range-----

mincolor=-1.0
maxcolor=1.0
;mincolor=min(log10rho)
;maxcolor=max(log10rho)


;-----set the number of colorbar level----

ncolors=64
colorlevels=mincolor+findgen(ncolors)*(maxcolor-mincolor)/ncolors
Hcolorbar=[[colorlevels], [colorlevels]]
Vcolorbar=transpose(Hcolorbar)


;-----set rho in range-----

for i=0, nx-1 do begin
	for j=0, ny-1 do begin

;-----spiral potential calcuration-----

spiralpot(j, i)=-(r(i)/7.0)^2/(1.0+(r(i)/7.0)^2)^(1.5)*cos(2*(phi(j)+alog(r(i))/tan(!pi/12.0)))

		if spiralpot(j, i) ge maxcolor then spiralpot(j, i)=colorlevels(ncolors-1)
		if spiralpot(j, i) le colorlevels(2) then spiralpot(j, i)=colorlevels(2)
endfor
endfor


;set_plot, 'ps'
;device, xsize=33, ysize=30, file='./densityeps/zplanerho-'+idlnumstring+'.eps', /color, /encapsulated

polar_contour, spiralpot, phi, r, nlevels=ncolors, xrange=[-15.0, 15.0], yrange=[-15.0, 15.0], zrange=[mincolor, maxcolor], position=[0.18, 0.18, 0.98, 0.98], xstyle=1, xtitle='x, [kpc]', ystyle=1, ytitle='y, [kpc]', /fill, charsize=2.5


;omega_sp=12.2km/s/kpc, i_sp=15*!pi/180
;for phasenum=-51, 51 do begin
;	oplot, exp(tan(!pi/12)*((phasenum)/2*!pi+phi+12.2e+5/kpc*1e+7*yr*idlnum)), phi, thick=3, color=0, /polar
;	phasenum=phasenum+1
;endfor

for phasenum=-10, 10 do begin
oplot, exp(tan(!pi/12.0)*(phasenum/2.0*!pi-phi)), phi, /polar

phasenum=phasenum+1
endfor

;contour,Vcolorbar, [0.0, 0.2], colorlevels, levels=colorlevels, position=[0.075, 0.18, 0.1, 0.48], xtickformat='(A1)', xstyle=7, xticks=0, xminor=0, ystyle=1, yticks=2, yminor=1, ytickv=[-26, -25, -24], charsize=3.0, /fill, /noerase

contour,Vcolorbar, [0.0, 0.2], colorlevels, levels=colorlevels, position=[0.075, 0.18, 0.1, 0.48], xtickformat='(A1)', xstyle=7, xticks=0, xminor=0, ystyle=1, yticks=4, yminor=1, ytickv=[-1.0, -0.5, 0, 0.5, 1.0], charsize=1.5, /fill, /noerase


arrow, 0.0, 0.0, 0.5, 0.25 ,/data


print, filename

;contour, z, x, y, levels=[0.5, 0.7, 0.9], position=[0.2, 0.1, 0.9, 0.9], /noerase, charsize=2, /follow


;write_png, './densityeps/xyplane-density'+idlnumstring+'.png', tvrd(true=1)


;device, /close


free_lun, lun


end


