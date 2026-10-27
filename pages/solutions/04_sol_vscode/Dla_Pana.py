import math


g = 9.81
v0 = 37.0
angle_deg = 52.0
angle = math.radians(angle_deg)

vx0 = v0 * math.cos(angle)
vy0 = v0 * math.sin(angle)
flight_time = 2.0 * vy0 / g
max_height = vy0**2 / (2.0 * g)
range_m = vx0 * flight_time

print(f"Predkosc poczatkowa: {v0:.2f} m/s")
print(f"Kat wyrzutu: {angle_deg:.1f} stopni")
print(f"Czas lotu: {flight_time:.2f} s")
print(f"Maksymalna wysokosc: {max_height:.2f} m")
print(f"Zasieg: {range_m:.2f} m")
print("\nCzas [s] | x [m] | y [m]")

for step in range(7):
	t = flight_time * step / 6
	x = vx0 * t
	y = vy0 * t - 0.5 * g * t**2
	print(f"{t:8.2f} | {x:6.2f} | {y:6.2f}")
