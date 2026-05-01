import math

Earth_mu = 3.986004418e5  # km^3/s^2
Earth_radius = 6378  # km
Altitude = 400  # km

def orbital_period(alt_km):
	a = 6378 + alt_km  # altitude du satellite
	T = 2*math.pi*math.sqrt(a**3/Earth_mu)   # calcul de la période orbitale

	return T/60

print(f"The orbital period at {Altitude} km is equal to {orbital_period(Altitude):.2f} minutes.")
