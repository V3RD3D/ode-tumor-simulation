# Proyecto EDO — Simulación numérica de modelos tumorales

![Python](https://img.shields.io/badge/Python-3.x-blue)
![NumPy](https://img.shields.io/badge/NumPy-2.x-blueviolet)
![Matplotlib](https://img.shields.io/badge/Matplotlib-3.x-orange)
![GitHub](https://img.shields.io/badge/GitHub-V3RD3D%2Fproyecto-181717?logo=github)

> Trabajo final de la materia **Ecuaciones Diferenciales Ordinarias**.  
> Basado en el artículo *Simple ODE Models of Tumor Growth and Anti-Angiogenic or Radiation Treatment* de **R. K. Sachs, L. R. Hlatky y P. Hahnfeldt** (2001).

---

##  Descripción

Este repositorio contiene la parte de **simulación numérica** del proyecto. Se integran dos sistemas de EDO con el método de **Euler explícito**, tal como pide la orientación de la materia:

1. **Angiogénesis con tratamiento antiangiogénico**  
   Reproducción cualitativa de la Figura 1 del artículo: el tumor crece, estimula la neovascularización, el fármaco reduce la capacidad de carga y el tumor se contrae durante el tratamiento.

2. **Modelo lineal-cuadrático (LQ) de radiobiología**  
   Simulación de una dosis aguda (`G = 1`) y comparación con la fórmula cerrada:

   ```text
   ln(N/N0) = -α·D - β·D²
   ```

Además se incluye una **traducción no oficial al español** del artículo original.

---

## Modelos simulados

### 1. Angiogénesis — ecuaciones (6)-(7)

```text
dV/dt = λ · V · ln(K/V)
dK/dt = -u(t)·K + ω·V - γ·V^(2/3)·K
```

donde:

- `V(t)` es el volumen tumoral.
- `K(t)` es la capacidad de carga dinámica, proporcional a la neovascularización.
- `u(t) = σ·c(t)` es la intensidad efectiva del fármaco.
- El esquema de tratamiento son bolos diarios en los días **10 a 19**, con decaimiento exponencial:

```text
u(t) = Σ  u0 · exp(-r·(t - t0)),   t0 = 10, 11, ..., 19
```

El equilibrio sin tratamiento satisface:

```text
V* = K* = (ω/γ)^(3/2) = K_inf
```

---

### 2. Radiobiología LQ — ecuaciones (8)-(9)

```text
dU/dt = δ·R(t) - U/τ - γ2·U²
dN/dt = -(α·R(t) + γ2·U²)·N
```

donde:

- `U(t)` es el número promedio de roturas de doble cadena (DSB) por célula.
- `N(t)` es el número de células clonogénicas.
- `R(t) = dD/dt` es la tasa de dosis.

Para una **dosis aguda** se aplica el salto exacto en `t = 0`:

```text
U(0+)  = δ·D
N(0+)  = N0 · exp(-α·D)
```

y luego se integra con `R(t) = 0`. El parámetro `γ2` se elige para recuperar el modelo LQ:

```text
β = δ² · γ2 · τ / 2
```

---

## ⚙️ Instalación y uso

Requiere Python 3 y las librerías `numpy` y `matplotlib`:

```powershell
pip install numpy matplotlib
```

Para ejecutar, primero entra a la carpeta donde clonaste el repositorio:

```powershell
git clone https://github.com/V3RD3D/proyecto.git
cd proyecto
pip install numpy matplotlib
python simulaciones.py
```

> Si ya renombraste el repositorio, reemplaza `proyecto` por el nombre nuevo en el `clone` y en el `cd`.

Si ya descargaste el proyecto, solo ubícate en la carpeta que contiene `simulaciones.py` y ejecuta:

```powershell
python simulaciones.py
```

En Windows también funciona:

```powershell
py simulaciones.py
```

El programa:

1. Imprime tablas en la consola.
2. Crea la carpeta `figuras/` si no existe.
3. Guarda dos gráficos en PNG a 200 dpi.

---

## 📊 Resultados

### Angiogénesis

| Día | `V(t)` | `K(t)` |
|---:|---:|---:|
| 0  | 0.1000 | 0.1000 |
| 10 | 8.6300 | 10.6637 |
| 15 | 5.8933 | 6.7884 |
| 19 | 5.0606 | 6.1980 |
| 25 | 9.7816 | 11.1578 |
| 30 | 11.5216 | 11.8287 |

### Supervivencia LQ (`G = 1`)

| Dosis (Gy) | EDO `ln(N/N0)` | LQ `ln(N/N0)` | Error |
|---:|---:|---:|---:|
| 0  | 0.000000  | 0.000000  | 0.000000 |
| 2  | -0.719880 | -0.720000 | +0.000120 |
| 5  | -2.247018 | -2.250000 | +0.002982 |
| 10 | -5.973307 | -6.000000 | +0.026693 |

La diferencia crece con la dosis porque la deducción del modelo LQ descarta el término no lineal en la ecuación de reparación. El artículo señala que la aproximación es buena para dosis de hasta aproximadamente **5 Gy**.

---

##  Figuras generadas

Las figuras se guardan automáticamente en `figuras/`:

- `fig1_angiogenesis.png` — tumor, capacidad de carga e intensidad del fármaco.
- `fig2_lq.png` — curva de supervivencia y reparación tras 2 Gy.

---

##  Parámetros principales

| Parámetro | Valor | Unidad | Comentario |
|---|---:|---|---|
| `lam` | 0.50 | 1/día | Tasa de crecimiento de Gompertz |
| `omega` | 5.85 | 1/día | Estimulación angiogénica |
| `K_inf` | 12.0 | u.a. | Capacidad de carga sin tratamiento |
| `gamma` | 1.1161 | 1/día | Calculado como `omega / K_inf^(2/3)` |
| `u0` | 6.0 | 1/día | Pico del antiangiogénico |
| `r_dec` | 2.0 | 1/día | Decaimiento del bolo |
| `V0 = K0` | 0.10 | u.a. | Condiciones iniciales |
| `T_ang` | 30 | días | Tiempo total simulado |
| `h_ang` | 0.001 | días | Paso de Euler |
| `alpha` | 0.30 | Gy⁻¹ | Daño lineal por radiación |
| `beta` | 0.03 | Gy⁻² | Daño cuadrático |
| `tau` | 0.50 | h | Tiempo de reparación |
| `delta` | 40.0 | DSB/Gy | Roturas dobles por unidad de dosis |
| `gamma2` | 7.5e-5 | 1/(DSB·h) | Calculado como `2·beta/(delta²·tau)` |

Los valores de `alpha` y `beta` son típicos en radiobiología (`alpha/beta = 10 Gy`); `delta` y `tau` se justifican en el artículo. Los parámetros de angiogénesis son ilustrativos, elegidos para reproducir cualitativamente la Figura 1.

---

##  Estructura del repositorio

```text
.
├── simulaciones.py
├── figuras/
│   ├── fig1_angiogenesis.png
│   └── fig2_lq.png
├── 5. Models of Anti-Angiogenic Treatment.pdf
├── 5. Modelos de Tratamiento Antiangiogenico (traduccion).pdf
├── 5. Modelos de Tratamiento Antiangiogenico (traduccion).tex
├── figura1_traduccion.png
├── Orientación Proyecto EDO.pdf
└── README.md
```

---

##  Traducción al español

El archivo `5. Modelos de Tratamiento Antiangiogenico (traduccion).pdf` contiene una traducción no oficial del artículo, con:

- resumen y palabras clave;
- secciones 1 a 6;
- ecuaciones (1)-(12) con la numeración original;
- Figura 1;
- las 57 referencias.

El original es un escaneo con capa OCR defectuosa, por lo que las ecuaciones fueron reconstruidas y dejadas consistentes con el modelo.

---

##  Topics sugeridos

```text
ordinary-differential-equations
numerical-simulation
numerical-methods
mathematical-biology
biomathematics
computational-biology
cancer-modeling
tumor-growth
tumor-dynamics
angiogenesis
anti-angiogenic-therapy
radiobiology
radiotherapy
linear-quadratic-model
euler-method
python
numpy
matplotlib
scientific-computing
open-science
```

---

## Créditos

- Artículo base: **R. K. Sachs, L. R. Hlatky & P. Hahnfeldt**, *Mathematical and Computer Modelling* **33** (2001), 1297-1305.
- Orientación del proyecto: documento `Orientación Proyecto EDO.pdf`.

Este es un trabajo académico. Si deseas reutilizarlo, menciona la fuente original.
