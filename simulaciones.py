# -*- coding: utf-8 -*-
"""
Proyecto final de EDO -- Articulo 5:
"Simple ODE Models of Tumor Growth and Anti-Angiogenic or Radiation Treatment"
(Sachs, Hlatky & Hahnfeldt, Math. Comput. Modelling 33 (2001) 1297-1305)

Parte 4 de la orientacion: SIMULACION NUMERICA
  1) Sistema de angiogenesis (ecs. 6-7) con tratamiento antiangiogenico
     -> reproduce cualitativamente la Fig. 1 del articulo.
  2) Sistema LQ de radiacion (ecs. 8-9) con G = 1 (dosis aguda)
     -> comparar con la formula cerrada  ln(N/N0) = -alpha*D - beta*D^2.

Metodo: Euler explicito (como pide la orientacion).
"""

import os
import numpy as np
import matplotlib.pyplot as plt



# 1) ANGIOGENESIS:  dV/dt = lambda*V*ln(K/V)
#                   dK/dt = -u(t)*K + omega*V - gamma*V^(2/3)*K
#    (ecs. (6)-(7); u(t) = sigma*c(t) es la intensidad efectiva del farmaco)


# --- Parametros ilustrativos (ajustables)
lam = 0.50     # 1/dia   tasa de crecimiento del tumor (Gompertz)
omega = 5.85     # 1/dia   estimulacion angiogenica producida por el tumor
K_inf = 12.0     # u.a.    capacidad de carga sin tratamiento = (omega/gamma)^(3/2)
gamma = omega / K_inf**(2.0/3.0)   # 1/dia  inhibicion de largo alcance
u0 = 6.0      # 1/dia   pico de u(t) = sigma*c0 (intensidad del antiangiogenico)
r_dec = 2.0      # 1/dia   decaimiento del bolo
V0 = 0.10     # volumen inicial del tumor
K0 = 0.10     # capacidad de carga inicial
t_inicio, t_fin = 10, 19            # tratamiento diario en los dias 10..19
T_ang = 30.0     # dias simulados
h_ang = 1e-3     # paso de Euler (dias)


def dosis_antiangiogenica(t):
    """u(t) = sigma*c(t), con un bolo diario c ~ c0*exp(-r*(t-t0))."""
    u = 0.0
    for t0 in range(t_inicio, t_fin + 1):
        if t >= t0:
            u += u0 * np.exp(-r_dec * (t - t0))
    return u


def simular_angiogenesis():
    n = int(T_ang / h_ang)
    t = np.linspace(0.0, T_ang, n + 1)
    V = np.zeros(n + 1)
    K = np.zeros(n + 1)
    u = np.zeros(n + 1)
    V[0], K[0] = V0, K0

    for i in range(n):
        u[i] = dosis_antiangiogenica(t[i])
        dV = lam * V[i] * np.log(max(K[i], 1e-12) / max(V[i], 1e-12))
        dK = -u[i] * K[i] + omega * V[i] - gamma * V[i]**(2.0/3.0) * K[i]
        V[i + 1] = max(V[i] + h_ang * dV, 1e-12)
        K[i + 1] = max(K[i] + h_ang * dK, 1e-12)
    u[-1] = dosis_antiangiogenica(t[-1])
    return t, V, K, u


def figura_angiogenesis(t, V, K, u, carpeta):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    ax1.axvspan(t_inicio, t_fin + 1, color="tab:red", alpha=0.10,
                label="tratamiento (dias 10-19)")
    ax1.plot(t, V, color="tab:blue", lw=2.2, label=r"tumor $V(t)$")
    ax1.plot(t, K, color="tab:green", lw=2.0, ls="--",
             label=r"capacidad de carga $K(t)$")
    ax1.set_ylabel("volumen (u.a.)")
    ax1.set_title("Angiogenesis + tratamiento antiangiogenico (ecs. 6-7)")
    ax1.legend(loc="lower right")
    ax1.grid(alpha=0.3)

    ax2.plot(t, u, color="tab:red", lw=1.6)
    ax2.set_xlabel("tiempo (dias)")
    ax2.set_ylabel(r"$u(t)=\sigma c(t)$ (1/dia)")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    ruta = os.path.join(carpeta, "fig1_angiogenesis.png")
    fig.savefig(ruta, dpi=200)
    plt.close(fig)
    return ruta



# 2) MODELO LQ:  dU/dt = delta*R(t) - U/tau - gamma2*U^2
#                dN/dt = -(alpha*R(t) + gamma2*U^2)*N
#    (ecs. (8)-(9));  G = 1  <=>  dosis aguda en t = 0


alpha  = 0.30      # Gy^-1
beta   = 0.03      # Gy^-2
tau    = 0.50      # h
delta  = 40.0      # DSB por Gy y por celula
gamma2 = 2.0 * beta / (delta**2 * tau)   # 1/(DSB*h);  beta = delta^2*gamma2*tau/2


def simular_dosis_aguda(D, T=24.0, h=1e-3):
    """Dosis aguda D (Gy) aplicada en t=0.

    La delta de Dirac se aplica como salto exacto:
        U(0+) = delta*D,      N(0+) = N0*exp(-alpha*D)
    y despues se integran (8)-(9) con R(t)=0 (reparacion/malreparacion).

    Devuelve ln(N(T)/N0).
    """
    U = delta * D
    logN = -alpha * D
    n = int(T / h)
    for _ in range(n):
        dU = -U / tau - gamma2 * U * U
        dlogN = -gamma2 * U * U
        U += h * dU
        logN += h * dlogN
    return logN


def figura_lq(carpeta):
    D = np.linspace(0.0, 10.0, 41)
    logN_ode = np.array([simular_dosis_aguda(d) for d in D])
    logN_lq = -alpha * D - beta * D**2

    #Figura: curva de supervivencia 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    ax1.plot(D, logN_lq, color="tab:red", lw=2.0,
             label=r"LQ: $-\alpha D-\beta D^2$")
    ax1.plot(D, logN_ode, "o", ms=3.5, color="tab:blue",
             label="EDOs (8)-(9), Euler")
    ax1.set_xlabel("dosis D (Gy)")
    ax1.set_ylabel(r"$\ln(N/N_0)$")
    ax1.set_title("Supervivencia celular (G = 1)")
    ax1.legend()
    ax1.grid(alpha=0.3)

    #Figura: evolucion temporal para D = 2 Gy 
    D0 = 2.0
    T = 8.0
    h = 1e-3
    n = int(T / h)
    tt = np.linspace(0.0, T, n + 1)
    U = np.zeros(n + 1)
    logN = np.zeros(n + 1)
    U[0] = delta * D0
    logN[0] = -alpha * D0
    for i in range(n):
        dU = -U[i] / tau - gamma2 * U[i] * U[i]
        dlogN = -gamma2 * U[i] * U[i]
        U[i + 1] = U[i] + h * dU
        logN[i + 1] = logN[i] + h * dlogN

    ax2.plot(tt, U, color="tab:green", lw=2.0, label=r"$U(t)$ (DSB/celula)")
    ax2.set_xlabel("tiempo (h)")
    ax2.set_ylabel("DSB por celula", color="tab:green")
    ax2.tick_params(axis="y", labelcolor="tab:green")
    ax2.grid(alpha=0.3)

    ax2b = ax2.twinx()
    ax2b.plot(tt, np.exp(logN), color="tab:blue", lw=2.0, ls="--",
              label=r"$N(t)/N_0$")
    ax2b.set_ylabel(r"$N(t)/N_0$", color="tab:blue")
    ax2b.tick_params(axis="y", labelcolor="tab:blue")

    ax2.set_title("Reparacion tras dosis aguda D = 2 Gy")

    fig.tight_layout()
    ruta = os.path.join(carpeta, "fig2_lq.png")
    fig.savefig(ruta, dpi=200)
    plt.close(fig)
    return D, logN_ode, logN_lq, ruta


# 
def main():
    carpeta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras")
    os.makedirs(carpeta, exist_ok=True)

    print("=" * 72)
    print("1) ANGIOGENESIS (ecs. 6-7)")
    print("=" * 72)
    t, V, K, u = simular_angiogenesis()
    ruta1 = figura_angiogenesis(t, V, K, u, carpeta)
    print(f"gamma = omega/K_inf^(2/3) = {gamma:.5f} 1/dia")
    print(f"equilibrio sin tratamiento (omega/gamma)^(3/2) = {K_inf:.3f}")
    for dia in (0, 10, 15, 19, 20, 25, 30):
        i = int(dia / h_ang)
        print(f"  t = {dia:5.1f} d   V = {V[i]:8.4f}   K = {K[i]:8.4f}   u = {u[i]:7.4f}")
    print("figura:", ruta1)

    print()
    print("=" * 72)
    print("2) MODELO LQ (ecs. 8-9), G = 1")
    print("=" * 72)
    print(f"gamma2 = 2*beta/(delta^2*tau) = {gamma2:.6e} 1/(DSB*h)")
    print(f"beta teorico = delta^2*gamma2*tau/2 = {delta**2*gamma2*tau/2:.6f} Gy^-2")
    D, logN_ode, logN_lq, ruta2 = figura_lq(carpeta)
    print("  D (Gy)    ln(N/N0) ODE    ln(N/N0) LQ     error")
    for k in range(0, len(D), 5):
        print(f"  {D[k]:6.2f}    {logN_ode[k]:12.6f}    {logN_lq[k]:12.6f}"
              f"{logN_ode[k]-logN_lq[k]:+.2e}")
    print("figura:", ruta2)
    print()
    print("Listo.")


if __name__ == "__main__":
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 12,
        "figure.autolayout": False,
    })
    main()
