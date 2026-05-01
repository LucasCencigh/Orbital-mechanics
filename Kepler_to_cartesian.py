import numpy as np 
import math

Earth_mu = 3.986004418e5  # km^3/s^2

def kepler_to_cartesian(a, e, i, omega, Omega, nu, mu=Earth_mu):
    ''' a : demi-grand axe (km)
        e : excentricité
        i : inclinaison, angle de nutation ou angle entre le plan orbital et le plan de référence (radians)
        omega : argument du périastre (radians)
        Omega : longitude du noeud ascendant (radians)
        nu : anomalie vraie (radians)
        mu : paramètre gravitationnel de la Terre (G*M_T) (km^3/s^2)
    '''
    # 1. Calcul du rayon de l'orbite r_norm et des vitesses radiale v_r et transversale v_t

    r_norm = (a*(1-e**2))/(1+e*np.cos(nu))

    v_r = np.sqrt(mu/(a*(1-e**2)))*e*np.sin(nu)

    v_t = np.sqrt(mu/(a*(1-e**2)))*(1+e*np.cos(nu))

    # 2. Calcul du rayon et de la vitesse dans le plan orbital

    r_orbital = np.array([r_norm*np.cos(nu), 
                          r_norm*np.sin(nu),
                          0])
    
    v_orbital = np.array([v_r*np.cos(nu)-v_t*np.sin(nu),
                          v_r*np.sin(nu)+v_t*np.cos(nu),
                          0])
    
    # 3. Calcul de la matrice de rotation R = Rz_omega * Rx_i * Rz_Omega

    #Matrice de rotation d'omega autour de z 
    Rz_omega = np.array([
        [np.cos(omega), -np.sin(omega), 0],
        [np.sin(omega), np.cos(omega), 0],
        [0,0,1]
    ])

    #Matrice de rotation de i autour de x
    Rx_i = np.array([
        [1,0,0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i), np.cos(i)],
    ])

    #Matrice de rotation de Omega autour de z
    Rz_Omega = np.array([
        [np.cos(Omega), -np.sin(Omega), 0],
        [np.sin(Omega), np.cos(Omega), 0],
        [0,0,1]
    ])

    #Calcul du double produit matriciel
    R = Rz_omega @ Rx_i @ Rz_Omega

    # 4. Calcul des vecteurs position et vitesse dans le référentiel inertiel (ECI : Earth-centered inertial)

    r_ECI = R @ r_orbital

    v_ECI = R @ v_orbital

    
    return r_ECI, v_ECI

a = 42164
e = 0
i = 0
omega = 0
Omega = 0
nu = 0

r, v = kepler_to_cartesian(a, e, i, Omega, omega, nu)

print(f"Position:{r}km")
print(f"Velocity magnitude : {v} km/s")



