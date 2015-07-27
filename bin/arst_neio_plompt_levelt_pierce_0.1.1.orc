sr = 96000
ksmps = 1
nchnls = 2

instr 1
	
ifrequency = p4
idb = p5
iamplitude = ampdb (idb)

; env table
ienvattack = p3*0.09
ienvdecay = p3*0.5
ienvlevel = 0.55
ienvrelease = p3*0.2
aenvelope mxadsr ienvattack , ienvdecay , ienvlevel , ienvrelease

asignal_l oscil iamplitude, ifrequency, 1; loscil iamplitude, ifrequency, 1, 52.63, 1
asignal_l = asignal_l * aenvelope

asignal_r oscil iamplitude, ifrequency, 2; loscil iamplitude, ifrequency, 2, 52.63, 1
asignal_r = asignal_r * aenvelope


outs asignal_l , asignal_r
endin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
instr 2

ifrequency = p4
idb = p5
iamplitude = ampdb (idb)

; env
;ienvattack = 0.05
;ienvdecay = 0.3
;ienvlevel = 0.4
;ienvrelease = 0.2
;aenvelope mxadsr ienvattack , ienvdecay , ienvlevel , ienvrelease

asignal_l loscil iamplitude, ifrequency, 1, 52.63, 1
;asignal_l = asignal_l * aenvelope

asignal_r loscil iamplitude, ifrequency, 2, 52.63, 1
;asignal_r = asignal_r * aenvelope

; damping
iattack = 0.01
irelease = 0.1
isustain = p3
p3 = iattack + isustain + irelease
adamping linsegr 0.0 , iattack , 1.0 , isustain , 1.0 , irelease , 0.0

asignal_l = asignal_l * adamping
asignal_r = asignal_r * adamping

outs asignal_l , asignal_r
endin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
instr 3

ifrequency = p4
idb = p5
iamplitude = ampdb (idb)

; env
;ienvattack = 0.05
;ienvdecay = 0.3
;ienvlevel = 0.4
;ienvrelease = 0.2
;aenvelope mxadsr ienvattack , ienvdecay , ienvlevel , ienvrelease

asignal_l loscil iamplitude, ifrequency, 3, 52.63
;asignal_l = asignal_l * aenvelope

asignal_r loscil iamplitude, ifrequency, 3, 52.63
;asignal_r = asignal_r * aenvelope

; damping
;iattack = 0.01
;irelease = 0.1
;isustain = p3
;p3 = iattack + isustain + irelease
;adamping linsegr 0.0 , iattack , 1.0 , isustain , 1.0 , irelease , 0.0

;asignal_l = asignal_l * adamping
;asignal_r = asignal_r * adamping

outs asignal_l , asignal_r
endin
