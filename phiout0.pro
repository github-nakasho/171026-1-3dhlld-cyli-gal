
pro phiout0


;-----set up window size, window title-----

window, 0, xsize=800, ysize=800, title='phi-rho'


;-----set device & colors-----

device, decompose=0
device, true_color=24
;!p.color=0
;!p.background=16777215
;!p.background=255
!p.multi=[0, 1, 4]

;---number of grid----

ny=64


;-----r, phi, rho array definition-----

phi=fltarr(ny)
rho=fltarr(ny)
vperp=fltarr(ny)
vpara=fltarr(ny)
pot=fltarr(ny)

;-----reading coordinates-----


openr, lun, './azimuth-0/phi.dat', /get_lun
readf, lun, phi


;-----start main loop-----

idlnum=0
for idlnum=0, 400 do begin


;-----change filename to output numbering-----

idlnumstring=string(idlnum, format='(i0)')
filename='./azimuth-0/potential-'+idlnumstring+'.dat'


;-----reading pressure-----

openr, lun, filename, /get_lun
readf, lun, pot


;-----change filename to output numbering-----

idlnumstring=string(idlnum, format='(i0)')
filename='./azimuth-0/rho-'+idlnumstring+'.dat'


;-----reading rho-----

openr, lun, filename, /get_lun
readf, lun, rho


;-----change filename to output numbering-----

idlnumstring=string(idlnum, format='(i0)')
filename='./azimuth-0/vperp-'+idlnumstring+'.dat'


;-----reading vperp-----

openr, lun, filename, /get_lun
readf, lun, vperp


;-----change filename to output numbering-----

idlnumstring=string(idlnum, format='(i0)')
filename='./azimuth-0/vpara-'+idlnumstring+'.dat'


;-----reading vpara-----

openr, lun, filename, /get_lun
readf, lun, vpara

;-----normarization-----

rho=rho*1.67e-24
vperp=vperp*6.558e+1
vpara=vpara*6.558e+1


;-----set logscale-----

log10rho=alog10(rho)


;-----set numbering-----

if idlnum lt 10 then begin
	idlnumstring='00'+idlnumstring
endif else begin
	if idlnum lt 100 then idlnumstring='0'+idlnumstring
endelse


;set_plot, 'ps'
;device, xsize=33, ysize=30, file='./benergyeps/zplanerho-'+idlnumstring+'.eps', /color, /encapsulated

plot, phi, log10rho, xstyle=1, yrange=[-25.0, -23.0]
plot, phi, vperp, xstyle=1, yrange=[0, 80]
plot, phi, vpara, xstyle=1, yrange=[75, 175]
plot, phi, pot, xstyle=1, yrange=[-0.05, 0.05]


print, filename
close, /all

;contour, z, x, y, levels=[0.5, 0.7, 0.9], position=[0.2, 0.1, 0.9, 0.9], /noerase, charsize=2, /follow


;write_png, './densityeps/xyplane-density'+idlnumstring+'.png', tvrd(true=1)


;device, /close


wait, 0.5

free_lun, lun


endfor
!p.multi=0

end


