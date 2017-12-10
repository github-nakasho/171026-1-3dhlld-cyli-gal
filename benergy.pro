
pro benergy


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
benergy=fltarr(ny, nx)


;-----reading coordinates-----

openr, lun, './idl0/r.txt', /get_lun
readf, lun, r

r=r*10


openr, lun, './idl0/phi.txt', /get_lun
readf, lun, phi


;-----start main loop-----

idlnum=0
for idlnum=0, 250 do begin


;-----change filename to output numbering-----

idlnumstring=string(idlnum, format='(i0)')
filename='./idl'+idlnumstring+'/zplanebenergy-'+idlnumstring+'.dat'


;-----reading pressure-----

openr, lun, filename, /get_lun
readf, lun, benergy


;-----normarization-----

benergy=benergy*7.18e-11


;-----set logscale-----

log10be=alog10(benergy)


;-----choosing color table-----

;loadct, 0    
;loadct, 1
;loadct, 2    
;loadct, 3
;loadct, 4
loadct, 5
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
;loadct, 17
;loadct, 18 ;only 41
;loadct, 22
;loadct, 23
;loadct, 24
;loadct, 25
;loadct, 33
;loadct, 34
;loadct, 38
;loadct, 39


;-----set value range-----

mincolor=-16
maxcolor=-8
;mincolor=min(log10rho)
;maxcolor=max(log10rho)


;-----set the number of colorbar level----

ncolors=65
colorlevels=mincolor+findgen(ncolors)*(maxcolor-mincolor)/ncolors
Hcolorbar=[[colorlevels], [colorlevels]]
Vcolorbar=transpose(Hcolorbar)


;-----set pr in range-----

for i=0, nx-1 do begin
   for j=0, ny-1 do begin 
      if log10be(j, i) ge colorlevels(ncolors-1) then log10be(j, i)=colorlevels(ncolors-1)
      if log10be(j, i) le colorlevels(1) then log10be(j, i)=colorlevels(2)
   endfor
endfor


;-----set numbering-----

if idlnum lt 10 then begin
	idlnumstring='00'+idlnumstring
endif else begin
	if idlnum lt 100 then idlnumstring='0'+idlnumstring
endelse


;set_plot, 'ps'
;device, xsize=31, ysize=30, file='./benergyeps/zplanerho-'+idlnumstring+'.eps', /color, /encapsulated

polar_contour, log10be, phi, r, nlevels=ncolors, xrange=[-20.0, 20.0], yrange=[-20.0, 20.0], zrange=[mincolor, maxcolor], position=[0.2, 0.1, 0.98, 0.88], xstyle=1, xtitle='x, [kpc]', ystyle=1, ytitle='y, [kpc]', /fill, charsize=2.5


;omega_sp=12.2km/s/kpc, i_sp=15*!pi/180
;for phasenum=-51, 51 do begin
;	oplot, exp(tan(!pi/12)*((phasenum)/2*!pi+phi+12.2e+5/kpc*1e+7*yr*idlnum)), phi, thick=3, color=0, /polar
;	phasenum=phasenum+1
;endfor


;contour,Vcolorbar, [0.0, 0.2], colorlevels, levels=colorlevels, position=[0.075, 0.18, 0.1, 0.48], xtickformat='(A1)', xstyle=7, xticks=0, xminor=0, ystyle=1, yticks=2, yminor=1, ytickv=[-26, -25, -24], charsize=3.0, /fill, /noerase

contour,Vcolorbar, [0.0, 0.2], colorlevels, levels=colorlevels, position=[0.075, 0.1, 0.1, 0.4], xtickformat='(A1)', xstyle=7, xticks=0, xminor=0, ystyle=1, yticks=4, yminor=1, ytickv=[-16, -14, -12, -10, -8], charsize=2.5, /fill, /noerase


timestring=string(idlnum*0.05*1.492e+2, format='(i0)')

;arrow, 0.0, 0.0, 0.5, 0.25 ,/data
xyouts, 0.4, 0.92, 't='+timestring+', [Myr]', charsize=4, /normal 


print, filename

;contour, z, x, y, levels=[0.5, 0.7, 0.9], position=[0.2, 0.1, 0.9, 0.9], /noerase, charsize=2, /follow


;write_png, './densityeps/xyplane-density'+idlnumstring+'.png', tvrd(true=1)


;device, /close


free_lun, lun


endfor


end


