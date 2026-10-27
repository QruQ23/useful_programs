from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


V0 = 37.0
ANGLE_DEG = 52.0
G = 9.81
MASS = 0.145
RADIUS = 0.0366
AIR_DENSITY = 1.225
DRAG_COEFFICIENT = 0.47
AREA = np.pi * RADIUS**2
DT = 0.005
OUTPUT_DIR = Path(__file__).parent / "files" / "rzut_ukosny"


def acceleration(state):
    _, _, vx, vy = state
    speed = np.hypot(vx, vy)
    drag_factor = 0.5 * AIR_DENSITY * DRAG_COEFFICIENT * AREA / MASS * speed
    return np.array([vx, vy, -drag_factor * vx, -G - drag_factor * vy])


def rk4_step(state, dt):
    k1 = acceleration(state)
    k2 = acceleration(state + 0.5 * dt * k1)
    k3 = acceleration(state + 0.5 * dt * k2)
    k4 = acceleration(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_with_drag(angle_deg, dt=DT):
    angle = np.radians(angle_deg)
    state = np.array([0.0, 0.0, V0 * np.cos(angle), V0 * np.sin(angle)])
    times = [0.0]
    states = [state.copy()]

    while len(times) < 200_000:
        next_state = rk4_step(state, dt)
        next_time = times[-1] + dt

        if next_state[1] < 0:
            fraction = state[1] / (state[1] - next_state[1])
            next_state = state + fraction * (next_state - state)
            next_state[1] = 0.0
            next_time = times[-1] + fraction * dt
            states.append(next_state)
            times.append(next_time)
            break

        states.append(next_state)
        times.append(next_time)
        state = next_state
    else:
        raise RuntimeError("Nie znaleziono chwili lądowania w limicie symulacji.")

    return np.asarray(times), np.asarray(states)


def simulate_without_drag():
    angle = np.radians(ANGLE_DEG)
    flight_time = 2 * V0 * np.sin(angle) / G
    time = np.linspace(0.0, flight_time, 600)
    vx = np.full_like(time, V0 * np.cos(angle))
    vy = V0 * np.sin(angle) - G * time
    x = vx * time
    y = V0 * np.sin(angle) * time - 0.5 * G * time**2
    return time, np.column_stack((x, y, vx, vy))


def save_figures():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    angle = np.radians(ANGLE_DEG)
    t_ideal, ideal = simulate_without_drag()
    t_drag, drag = simulate_with_drag(ANGLE_DEG)

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ideal[:, 0], ideal[:, 1], label="Bez oporu")
    ax.plot(drag[:, 0], drag[:, 1], label="Z oporem powietrza")
    ax.set(title="Tor ruchu", xlabel="Położenie x [m]", ylabel="Wysokość y [m]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "tor_ruchu.svg")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].plot(t_ideal, ideal[:, 0], label="Bez oporu")
    axes[0].plot(t_drag, drag[:, 0], label="Z oporem")
    axes[0].set_ylabel("Położenie x [m]")
    axes[0].legend()
    axes[1].plot(t_ideal, ideal[:, 1], label="Bez oporu")
    axes[1].plot(t_drag, drag[:, 1], label="Z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="Wysokość y [m]")
    axes[1].legend()
    fig.suptitle("Położenie w funkcji czasu")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "polozenie_czas.svg")
    plt.close(fig)

    ideal_speed = np.hypot(ideal[:, 2], ideal[:, 3])
    drag_speed = np.hypot(drag[:, 2], drag[:, 3])
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].plot(t_ideal, ideal[:, 2], label="Bez oporu")
    axes[0].plot(t_drag, drag[:, 2], label="Z oporem")
    axes[0].set_ylabel("Prędkość vx [m/s]")
    axes[0].legend()
    axes[1].plot(t_ideal, ideal[:, 3], label="vy bez oporu")
    axes[1].plot(t_drag, drag[:, 3], label="vy z oporem")
    axes[1].plot(t_ideal, ideal_speed, "--", label="Szybkość bez oporu")
    axes[1].plot(t_drag, drag_speed, "--", label="Szybkość z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="Prędkość [m/s]")
    axes[1].legend(ncol=2, fontsize=8)
    fig.suptitle("Prędkość w funkcji czasu")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "predkosc_czas.svg")
    plt.close(fig)

    ideal_kinetic = 0.5 * MASS * (ideal[:, 2] ** 2 + ideal[:, 3] ** 2)
    ideal_potential = MASS * G * ideal[:, 1]
    drag_kinetic = 0.5 * MASS * (drag[:, 2] ** 2 + drag[:, 3] ** 2)
    drag_potential = MASS * G * drag[:, 1]
    drag_mechanical = drag_kinetic + drag_potential
    drag_power = 0.5 * AIR_DENSITY * DRAG_COEFFICIENT * AREA * drag_speed**3
    dissipated = np.concatenate(
        ([0.0], np.cumsum(0.5 * (drag_power[:-1] + drag_power[1:]) * np.diff(t_drag)))
    )

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(t_ideal, ideal_kinetic, label="Kinetyczna, bez oporu")
    axes[0].plot(t_ideal, ideal_potential, label="Potencjalna, bez oporu")
    axes[0].plot(t_ideal, ideal_kinetic + ideal_potential, "--", label="Mechaniczna, bez oporu")
    axes[0].set_ylabel("Energia [J]")
    axes[0].legend(fontsize=8)
    axes[1].plot(t_drag, drag_kinetic, label="Kinetyczna, z oporem")
    axes[1].plot(t_drag, drag_potential, label="Potencjalna, z oporem")
    axes[1].plot(t_drag, drag_mechanical, label="Mechaniczna, z oporem")
    axes[1].plot(t_drag, dissipated, "--", label="Rozproszona przez opór")
    axes[1].set(xlabel="Czas t [s]", ylabel="Energia [J]")
    axes[1].legend(fontsize=8)
    fig.suptitle("Bilans energii")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "energia.svg")
    plt.close(fig)

    angles = np.arange(10.0, 86.0, 5.0)
    ranges_ideal = V0**2 * np.sin(2 * np.radians(angles)) / G
    ranges_drag = np.array([simulate_with_drag(a)[1][-1, 0] for a in angles])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(angles, ranges_ideal, "o-", label="Bez oporu")
    ax.plot(angles, ranges_drag, "s-", label="Z oporem")
    ax.axvline(ANGLE_DEG, color="0.35", linestyle=":", label=f"Kąt początkowy {ANGLE_DEG:.0f}°")
    ax.set(title="Zasięg a kąt wyrzutu", xlabel="Kąt wyrzutu [°]", ylabel="Zasięg [m]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "zasieg_kat.svg")
    plt.close(fig)

    ideal_range = ideal[-1, 0]
    ideal_height = np.max(ideal[:, 1])
    drag_range = drag[-1, 0]
    drag_height = np.max(drag[:, 1])
    print(f"Bez oporu: czas lotu {t_ideal[-1]:.3f} s, zasięg {ideal_range:.3f} m, wysokość {ideal_height:.3f} m")
    print(f"Z oporem: czas lotu {t_drag[-1]:.3f} s, zasięg {drag_range:.3f} m, wysokość {drag_height:.3f} m")
    print(f"Wykresy zapisano w: {OUTPUT_DIR}")


if __name__ == "__main__":
    save_figures()