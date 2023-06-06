import numpy as np
import matplotlib.pyplot as plt
import math

#OMEGA = longitude do nodo ascendente
#incli vetor que armazena a inclinação de cada planeta em relação ao plano de orbita terreste (i)
#omg_coord vetor que armazena o argumento de periélio de cada planeta (w)
#semi_eix_m valor do semi eixo maior de rotação do planeta (a)
#exc_orb excentricidade orbital (e)
#an_med_plan anomalia média 
#an_verd_plan anomalia verdadeira
#dist_r distancia do planeta
#lamb_coord coordenada eliptica lambda
#beta_coord coordenada elíptica beta

'''            COMEMO DO PROGRAMA DE FATO               '''
## ano(y), mês(m), dia(d), hora (h) e minutos (m) dado a mim
y = 2117.
m = 7.
d = 24.
h = 2.
mim = 45.
## cálculo do tempo (t)
t_ym = 367*y - int((7*(y +(m + 9)/12))/4)+int((275*m)/9)+ d - 730530
t = t_ym + h/24 + mim/1440
t = int(t)
print(t)

###########################################################
''' Elementos Orbitais'''
omg_coord = np.zeros(8)
omg_coord[0] = 29.1241 + (1.01444*(10**-5)*t)
omg_coord[1] = 54.891 + (1.38374 * (10**-5 )*t)
omg_coord[2] = 282.9404 + (4.70935*(10**-5)*t)
omg_coord[3] = 286.5016 + (2.92961 *(10**-5)*t)
omg_coord[4] = 273.8777 + (1.64505 * (10**-5)*t)
omg_coord[5] = 339.3939 + (2.97661 *(10**-5)*t)
omg_coord[6] = 96.6612 + (3.0565 * (10**-5)*t)
omg_coord[7] = 272.8461 - (6.027 *(10**-5)*t)
 
incli = np.zeros(8)
incli[0] = (7.0047 + (5*(10**-8)*t))
incli[1] = (3.3946 + (2.75 * (10**-8) *t))
incli[2] = 0
incli[3] = (1.8497 -(1.78 * (10**-8)* t))
incli[4] = (1.303 - (1.557 * (10**-7) * t))
incli[5] = (2.4886 - (1.081 * (10**-7) *t))
incli[6] = (0.7733 + (1.9 *( 10**-8) *t))
incli[7] = (1.77 -(2.55 * (10**-7 )*t))

semi_eix_m = np.zeros(8)
semi_eix_m[0] = 0.387098
semi_eix_m[1] = 0.72333
semi_eix_m[2]  = 1
semi_eix_m[3]  = 1.523688
semi_eix_m[4]  = 5.20256
semi_eix_m[5]  = 9.55475
semi_eix_m[6]  = 19.18171 - 1.55*(10**-8)*t
semi_eix_m[7]  = 30.05826 + 3.313*(10**-8)*t

exc_orb = np.zeros(8)
exc_orb[0] = 0.205635 + (5.59 * (10**-10) * t)
exc_orb[1] = 0.006773 - 1.302*(10**-9)*t
exc_orb[2] = 0.016709 - 1.151*(10**-9)*t
exc_orb[3] = 0.093405 + 2.516*(10**-9)*t
exc_orb[4] = 0.048498 + 4.469*(10**-9)*t
exc_orb[5] = 0.055546 - 9.499*(10**-9)*t
exc_orb[6] = 0.047318 + 7.45*(10**-9)*t 
exc_orb[7] = 0.008606 + 2.15*(10**-9)*t 

an_med_plan = np.zeros(8)
an_med_plan[0] = (168.6562 + 4.0923344368*t)%360
an_med_plan[1] = (48.0052 + 1.6021302244*t)%360
an_med_plan[2] = (356.047 + 0.9856002585*t)%360
an_med_plan[3] = (18.6021 + 0.5240207766*t)%360
an_med_plan[4] = (19.895 + 0.0830853001*t)%360
an_med_plan[5] = (316.967 + 0.0334442282*t)%360
an_med_plan[6] = (142.5905 + 0.011725806*t)%360
an_med_plan[7] = (260.2471 + 0.005995147*t)%360 
############################################################ 
'''declaração dos arrays que receberão as coordenadas'''
an_verd_plan = np.zeros(8)
dist_r = np.zeros(8)
lamb_coord = np.zeros(8)
beta_coord = np.zeros(8)
intera = [0,1,2,3,4,5,6,7]
############################################################
erro = 5*10**-6 ##erro da anomalia média
############################################################
'''parâmetros do plot'''
np.random.seed(19680801)
img = plt.figure(figsize = (15,20))
ax = img.add_subplot(projection='polar', facecolor = 'black')
nomes = ['Marcúrio', 'Vênus', 'Terra', 'Marte', 'Jupter', 'Saturno', 'Netuno', 'Urano']
cores = ['peru','crimson','cornflowerblue','orangered','tan','yellow','slateblue','darkturquoise']
#############################################################

for i, nome, cor in zip(intera, nomes, cores):
    mean_E = math.radians(an_med_plan[i])
    delta_E = (math.radians(an_med_plan[i])- mean_E + exc_orb[i]*np.sin(mean_E))/(1. - exc_orb[i]*np.cos(mean_E))
    if delta_E/erro>1 or delta_E/erro<-1 :
        while delta_E>erro or delta_E<-erro :
            mean_E = mean_E + delta_E
            delta_E = (math.radians(an_med_plan[i]) - mean_E + exc_orb[i]*np.sin(mean_E))/(1. - exc_orb[i]*np.cos(mean_E))        
    mean_E = mean_E + delta_E    
    an_verd_plan[i] = np.arctan(2*np.sqrt((1+ math.radians(exc_orb[0])/(1-math.radians(exc_orb[0]))*np.tan(mean_E/2))))
    dist_r[i] = (semi_eix_m[i]*(1 - exc_orb[i]**2))/(1 + exc_orb[i]*np.cos(an_verd_plan[i]))
    lamb_coord[i] =  np.arctan(np.tan(math.radians(omg_coord[i])+an_verd_plan[i])*np.cos(math.radians(incli[i])))
    beta_coord[i] = np.arccos((np.cos(math.radians(omg_coord[i]) + an_verd_plan[i]))/(np.cos(lamb_coord[i])))
    if np.cos(math.radians(omg_coord[i]))<0:
        beta_coord[i] = beta_coord[i]+np.pi
    ax.scatter(math.degrees(lamb_coord[i]),dist_r[i], s=65, alpha=0.75, color = cor, label = nome)
    print(round(math.degrees(lamb_coord[i]), 2),round(math.degrees(beta_coord[i]),2))


angle = np.deg2rad(45)
ax.legend(loc='upper left',  bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
plt.show()
'''
img2 = plt.figure()
ax2 = img2.add_subplot(111, projection="hammer")
ax2.scatter(lamb_coord, beta_coord)
plt.grid(True)
plt.show'''