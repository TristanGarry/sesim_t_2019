######################
# microcosm food web #
######################

# species body sizes
B.subtilis   = 10
B.cereus     = 4
S.marcescens = 2
P.aurelia    = (111.6)
P.bursaria   = (101.3)
Spirostomum  = (843.8)
Colpidium    = (81.0)
T.pyriformis = (26.7)
Vorticella   = (79.7)
Monostyla    = (127.3)
Euplotes     = (85.4)
Blepharisma  = (471.3)

# interactions
# effect of collumn species (.c) on row species (r.)

# constants
o = 0.0007
p = 0.66
e = 0.8

# competition
bintra  = 0
spintra = 0
spinter = -0.01

# B. subtilis
i1.1  = bintra
i1.2  = spinter
i1.3  = spinter
i1.4  = -o*(P.aurelia   / B.subtilis)^p
i1.5  = -o*(P.bursaria  / B.subtilis)^p
i1.6  = -o*(Spirostomum / B.subtilis)^p
i1.7  = -o*(Colpidium   / B.subtilis)^p
i1.8  = -o*(T.pyriformis/ B.subtilis)^p
i1.9  = -o*(Vorticella  / B.subtilis)^p
i1.10 = -o*(Monostyla   / B.subtilis)^p
i1.11 = -o*(Euplotes    / B.subtilis)^p
i1.12 = -o*(Blepharisma / B.subtilis)^p

# B. cereus
i2.1  = spinter
i2.2  = bintra
i2.3  = spinter 
i2.4  = -o*(P.aurelia   / B.cereus)^p
i2.5  = -o*(P.bursaria  / B.cereus)^p
i2.6  = -o*(Spirostomum / B.cereus)^p
i2.7  = -o*(Colpidium   / B.cereus)^p
i2.8  = -o*(T.pyriformis/ B.cereus)^p
i2.9  = -o*(Vorticella  / B.cereus)^p
i2.10 = -o*(Monostyla   / B.cereus)^p
i2.11 = -o*(Euplotes    / B.cereus)^p
i2.12 = -o*(Blepharisma / B.cereus)^p

# S. marcescens
i3.1  = spinter
i3.2  = spinter
i3.3  = bintra
i3.4  = -o*(P.aurelia   / S.marcescens)^p
i3.5  = -o*(P.bursaria  / S.marcescens)^p
i3.6  = -o*(Spirostomum / S.marcescens)^p
i3.7  = -o*(Colpidium   / S.marcescens)^p
i3.8  = -o*(T.pyriformis/ S.marcescens)^p
i3.9  = -o*(Vorticella  / S.marcescens)^p
i3.10 = -o*(Monostyla   / S.marcescens)^p
i3.11 = -o*(Euplotes    / S.marcescens)^p
i3.12 = -o*(Blepharisma / S.marcescens)^p

# P. aurelia
i4.1  = abs(e*(-o*(P.aurelia / B.subtilis)^p))
i4.2  = abs(e*(-o*(P.aurelia / B.cereus)^p))
i4.3  = abs(e*(-o*(P.aurelia / S.marcescens)^p))
i4.4  = spintra
i4.5  = spinter
i4.6  = spinter
i4.7  = spinter
i4.8  = spinter
i4.9  = spinter
i4.10 = spinter
i4.11 = spinter
i4.12 = spinter

# P. bursaria
i5.1  = abs(e*(-o*(P.bursaria / B.subtilis)^p))
i5.2  = abs(e*(-o*(P.bursaria / B.cereus)^p))
i5.3  = abs(e*(-o*(P.bursaria / S.marcescens)^p))
i5.4  = spinter
i5.5  = spintra
i5.6  = spinter
i5.7  = spinter
i5.8  = spinter
i5.9  = spinter
i5.10 = spinter
i5.11 = spinter
i5.12 = spinter

# Spirostomum sp.
i6.1  = abs(e*(-o*(Spirostomum / B.subtilis)^p))
i6.2  = abs(e*(-o*(Spirostomum / B.cereus)^p))
i6.3  = abs(e*(-o*(Spirostomum / S.marcescens)^p))
i6.4  = spinter
i6.5  = spinter
i6.6  = spintra
i6.7  = spinter
i6.8  = spinter
i6.9  = spinter
i6.10 = spinter
i6.11 = spinter
i6.12 = spinter

# Colpidium sp.
i7.1  = abs(e*(-o*(Colpidium / B.subtilis)^p))
i7.2  = abs(e*(-o*(Colpidium / B.cereus)^p))
i7.3  = abs(e*(-o*(Colpidium / S.marcescens)^p))
i7.4  = spinter
i7.5  = spinter
i7.6  = spinter
i7.7  = spintra
i7.8  = spinter
i7.9  = spinter
i7.10 = spinter
i7.11 = -o*(Euplotes / Colpidium)^p
i7.12 = -o*(Blepharisma / Colpidium)^p

# Tetrahymena pyriformis
i8.1  = abs(e*(-o*(T.pyriformis / B.subtilis)^p))
i8.2  = abs(e*(-o*(T.pyriformis / B.cereus)^p))
i8.3  = abs(e*(-o*(T.pyriformis / S.marcescens)^p))
i8.4  = spinter
i8.5  = spinter
i8.6  = spinter
i8.7  = spinter
i8.8  = spintra
i8.9  = spinter
i8.10 = spinter
i8.11 = -o*(Euplotes / T.pyriformis)^p
i8.12 = -o*(Blepharisma / T.pyriformis)^p

# Vorticella sp.
i9.1  = abs(e*(-o*(Vorticella / B.subtilis)^p))
i9.2  = abs(e*(-o*(Vorticella / B.cereus)^p))
i9.3  = abs(e*(-o*(Vorticella / S.marcescens)^p))
i9.4  = spinter
i9.5  = spinter
i9.6  = spinter
i9.7  = spinter
i9.8  = spinter
i9.9  = spintra
i9.10 = spinter
i9.11 = -o*(Euplotes / Vorticella)^p
i9.12 = -o*(Blepharisma / Vorticella)^p

# Monostyla sp.
i10.1  = abs(e*(-o*(Monostyla / B.subtilis)^p))
i10.2  = abs(e*(-o*(Monostyla / B.cereus)^p))
i10.3  = abs(e*(-o*(Monostyla / S.marcescens)^p))
i10.4  = spinter
i10.5  = spinter
i10.6  = spinter
i10.7  = spinter
i10.8  = spinter
i10.9  = spinter
i10.10 = spintra
i10.11 = spinter
i10.12 = -o*(Blepharisma / Monostyla)^p

# Euplotes
i11.1  = abs(e*(-o*(Euplotes / B.subtilis)^p))
i11.2  = abs(e*(-o*(Euplotes / B.cereus)^p))
i11.3  = abs(e*(-o*(Euplotes / S.marcescens)^p))
i11.4  = spinter
i11.5  = spinter
i11.6  = spinter
i11.7  = abs(e*(-o*(Euplotes / Colpidium)^p))
i11.8  = abs(e*(-o*(Euplotes / T.pyriformis)^p))
i11.9  = abs(e*(-o*(Euplotes / Vorticella)^p))
i11.10 = spinter
i11.11 = spintra
i11.12 = -o*(Blepharisma / Euplotes)^p

# Blepharisma
i12.1  = abs(e*(-o*(Blepharisma / B.subtilis)^p))
i12.2  = abs(e*(-o*(Blepharisma / B.cereus)^p))
i12.3  = abs(e*(-o*(Blepharisma / S.marcescens)^p))
i12.4  = spinter  
i12.5  = spinter
i12.6  = spinter
i12.7  = abs(e*(-o*(Blepharisma / Colpidium)^p))
i12.8  = abs(e*(-o*(Blepharisma / T.pyriformis)^p))
i12.9  = abs(e*(-o*(Blepharisma / Vorticella)^p))
i12.10 = abs(e*(-o*(Blepharisma / Monostyla)^p))
i12.11 = abs(e*(-o*(Blepharisma / Euplotes)^p))
i12.12 = abs(e*(-o*(Blepharisma / Blepharisma)^p)) + (-o*(Blepharisma / Blepharisma)^p)

# assembling food-web
FW <- (rbind(cbind(i1.1,i1.2,i1.3,i1.4,i1.5,i1.6,i1.7,i1.8,i1.9,i1.10,i1.11,i1.12),
             cbind(i2.1,i2.2,i2.3,i2.4,i2.5,i2.6,i2.7,i2.8,i2.9,i2.10,i2.11,i2.12),
             cbind(i3.1,i3.2,i3.3,i3.4,i3.5,i3.6,i3.7,i3.8,i3.9,i3.10,i3.11,i3.12),
             cbind(i4.1,i4.2,i4.3,i4.4,i4.5,i4.6,i4.7,i4.8,i4.9,i4.10,i4.11,i4.12),
             cbind(i5.1,i5.2,i5.3,i5.4,i5.5,i5.6,i5.7,i5.8,i5.9,i5.10,i5.11,i5.12),
             cbind(i6.1,i6.2,i6.3,i6.4,i6.5,i6.6,i6.7,i6.8,i6.9,i6.10,i6.11,i6.12),
             cbind(i7.1,i7.2,i7.3,i7.4,i7.5,i7.6,i7.7,i7.8,i7.9,i7.10,i7.11,i7.12),
             cbind(i8.1,i8.2,i8.3,i8.4,i8.5,i8.6,i8.7,i8.8,i8.9,i8.10,i8.11,i8.12),
             cbind(i9.1,i9.2,i9.3,i9.4,i9.5,i9.6,i9.7,i9.8,i9.9,i9.10,i9.11,i9.12),
             cbind(i10.1,i10.2,i10.3,i10.4,i10.5,i10.6,i10.7,i10.8,i10.9,i10.10,i10.11,i10.12),
             cbind(i11.1,i11.2,i11.3,i11.4,i11.5,i11.6,i11.7,i11.8,i11.9,i11.10,i11.11,i11.12),
             cbind(i12.1,i12.2,i12.3,i12.4,i12.5,i12.6,i12.7,i12.8,i12.9,i12.10,i12.11,i12.12)))

#######
# end #
#######