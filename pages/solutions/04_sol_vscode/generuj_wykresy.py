"""Generuje wykresy do raportu o rzucie ukośnym.

Uruchom w katalogu tego pliku poleceniem: python generuj_wykresy.py
Wymagane pakiety: numpy i matplotlib.
"""

import matplotlib
from pathlib import Path

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).parent / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

V0 = 37.0
ALPHA = np.deg2rad(52.0)
G = 9.81
MASS = 0.145
RHO = 1.225
CD = 0.47
DIAMETER = 0.073
AREA = np.pi * (DIAMETER / 2) ** 2
DT = 0.002

BLUE = "#2463A6"
ORANGE = "#D47722"
GREEN = "#168578"
PURPLE = "#8455A4"
INK = "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def acceleration(state):
    vx, vy = state[2], state[3]
    speed = np.hypot(vx, vy)
    drag_per_mass = 0.5 * RHO * CD * AREA * speed / MASS
    return np.array([vx, vy, -drag_per_mass * vx,
                     -G - drag_per_mass * vy])


def integrate_with_drag():
    state = np.array([
        0.0, 0.0, V0 * np.cos(ALPHA), V0 * np.sin(ALPHA)
    ])
    states = [state.copy()]
    times = [0.0]
    time = 0.0

    while time < 30.0:
        k1 = acceleration(state)
        k2 = acceleration(state + 0.5 * DT * k1)
        k3 = acceleration(state + 0.5 * DT * k2)
        k4 = acceleration(state + DT * k3)
        next_state = state + DT * (k1 + 2*k2 + 2*k3 + k4) / 6
        next_time = time + DT

        if next_state[1] < 0:
            fraction = state[1] / (state[1] - next_state[1])
            landing_state = state + fraction * (next_state - state)
            states.append(landing_state)
            times.append(time + fraction * DT)
            break

        states.append(next_state.copy())
        times.append(next_time)
        state, time = next_state, next_time
    else:
        raise RuntimeError("Nie znaleziono chwili upadku w limicie 30 s.")

    return np.asarray(times), np.asarray(states)


def save_figure(fig, filename):
    fig.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


vx0 = V0 * np.cos(ALPHA)
vy0 = V0 * np.sin(ALPHA)
flight_time = 2 * vy0 / G
vacuum_range = vx0 * flight_time
vacuum_height = vy0**2 / (2 * G)
time_vacuum = np.linspace(0, flight_time, 1601)
vacuum = np.column_stack((
    vx0 * time_vacuum,
    vy0 * time_vacuum - 0.5 * G * time_vacuum**2,
    np.full_like(time_vacuum, vx0),
    vy0 - G * time_vacuum,
))

time_drag, drag = integrate_with_drag()
drag_range = drag[-1, 0]
drag_height = np.max(drag[:, 1])
drag_flight_time = time_drag[-1]

assert np.isclose(vx0**2 + vy0**2, V0**2)
assert abs(vacuum[-1, 1]) < 1e-10
assert drag_range < vacuum_range
assert drag_height < vacuum_height
assert drag_flight_time < flight_time

fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
ax.plot(vacuum[:, 0], vacuum[:, 1], color=BLUE, label="Bez oporu")
ax.plot(drag[:, 0], drag[:, 1], color=ORANGE, label="Z oporem $v^2$")
ax.scatter([vacuum_range, drag_range], [0, 0], color=[BLUE, ORANGE], zorder=3)
ax.scatter([vacuum[np.argmax(vacuum[:, 1]), 0],
             drag[np.argmax(drag[:, 1]), 0]],
           [vacuum_height, drag_height], color=[BLUE, ORANGE], zorder=3)
ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]",
       title="Tor lotu: model idealny i opór aerodynamiczny")
ax.legend()
save_figure(fig, "tor_lotu.png")

fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=False,
                         layout="constrained")
for values, label, color in [
    (vacuum, "Bez oporu", BLUE), (drag, "Z oporem", ORANGE)
]:
    times = time_vacuum if values is vacuum else time_drag
    axes[0].plot(times, values[:, 0], color=color, label=label)
    axes[1].plot(times, values[:, 1], color=color, label=label)
axes[0].set(ylabel="Położenie x [m]", title="Położenie w funkcji czasu")
axes[1].set(xlabel="Czas t [s]", ylabel="Wysokość y [m]")
for ax in axes:
    ax.legend()
save_figure(fig, "polozenie_czas.png")

fig, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True,
                         layout="constrained")
for values, times, label, color in [
    (vacuum, time_vacuum, "Bez oporu", BLUE),
    (drag, time_drag, "Z oporem", ORANGE),
]:
    vx, vy = values[:, 2], values[:, 3]
    axes[0].plot(times, vx, color=color, label=label)
    axes[1].plot(times, vy, color=color, label=label)
    axes[2].plot(times, np.hypot(vx, vy), color=color, label=label)
axes[0].set_ylabel("Prędkość $v_x$ [m/s]")
axes[1].set_ylabel("Prędkość $v_y$ [m/s]")
axes[2].set_ylabel("Szybkość $|v|$ [m/s]")
axes[2].set_xlabel("Czas t [s]")
axes[0].set_title("Składowe prędkości i szybkość")
for ax in axes:
    ax.legend(ncol=2)
save_figure(fig, "predkosc_czas.png")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5),
                         layout="constrained")
for ax, values, times, title, color in [
    (axes[0], vacuum, time_vacuum, "Bez oporu", BLUE),
    (axes[1], drag, time_drag, "Z oporem", ORANGE),
]:
    kinetic = 0.5 * MASS * (values[:, 2]**2 + values[:, 3]**2)
    potential = MASS * G * values[:, 1]
    ax.plot(times, kinetic, color=color, label="Kinetyczna")
    ax.plot(times, potential, color=GREEN, label="Potencjalna")
    ax.plot(times, kinetic + potential, color=PURPLE, ls="--",
            label="Mechaniczna")
    ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", title=title)
    ax.legend(fontsize=8)
fig.suptitle("Energia mechaniczna: zachowanie i straty")
save_figure(fig, "energia.png")

print(f"v0 = {V0:.2f} m/s; alfa = {np.rad2deg(ALPHA):.1f} deg")
print(f"Bez oporu: czas = {flight_time:.3f} s, zasięg = {vacuum_range:.3f} m, "
      f"wysokość maksymalna = {vacuum_height:.3f} m")
print(f"Z oporem: czas = {drag_flight_time:.3f} s, zasięg = {drag_range:.3f} m, "
      f"wysokość maksymalna = {drag_height:.3f} m")
print(f"Zapisano 4 wykresy w: {OUT}")