# Laboratorio de Factorización LU · Unidad 2

Solución de sistemas de ecuaciones lineales **A·x = b** por factorización: **Doolittle, Crout y Cholesky**,
en una app web y en Python.

Universidad Peruana Unión · Ingeniería de Sistemas · Fernando Istaña Flores
Curso: Métodos Numéricos · Docente: Jorge Luis Manrique Plasencia

---

## Qué hay aquí

| Archivo | Contenido |
|---|---|
| `index.html` | La app web. Un solo archivo, sin dependencias: se abre con doble clic o desde GitHub Pages |
| `lu_sel.py` | Los tres métodos en Python puro, más las sustituciones y el caso de la guía resuelto |
| `salida.txt` | La salida del programa, para pegarla en el informe |
| `informes/` | El informe de la GAA de la Sesión 6 en PDF |

## La app web

- **Editor del sistema** con n ajustable de 2 a 8. Cada ecuación que agregas suma también una incógnita,
  porque la factorización solo existe para matrices cuadradas.
- **Cómo se resuelve**: cada paso k muestra el multiplicador o la fórmula aplicada y la tabla resultante.
  En Doolittle es la matriz aumentada `[A | b]` con el pivote resaltado; en Crout y Cholesky, las matrices
  L y U llenándose columna por columna, con las entradas pendientes atenuadas.
- **Métodos**
  - *Doolittle* — L con diagonal 1, guarda los multiplicadores de la eliminación gaussiana. Admite pivoteo parcial (PA = LU).
  - *Crout* — la diagonal unitaria va en U; L se llena por columnas y U por filas, sin eliminación previa.
  - *Cholesky* — A = L·Lᵀ, solo para matrices simétricas y definidas positivas; cuesta la mitad que LU.
    La app avisa si A no es simétrica o si el radicando se vuelve negativo.
- **Varios vectores b**: agregas los que quieras y todos se resuelven con la misma factorización.
  El contador compara los flops de LU contra recalcular la eliminación en cada vector.
- **Verificación** en todo momento: det(A), ‖LU − A‖∞ y el residuo A·x − b por ecuación.

## Ejecutar la versión en Python

```bash
python lu_sel.py
```

Solo usa la biblioteca estándar (`math`). Las funciones son independientes:

```python
from lu_sel import doolittle, crout, cholesky, resolver

L, U = doolittle(A)          # una sola vez, O(n^3)
y, x = resolver(L, U, b)     # cada vector, O(n^2)
```

`resolver` sirve para los tres métodos porque la sustitución progresiva divide entre `L[i][i]`
(en Doolittle ese valor es 1, en Crout y Cholesky no).

## Caso resuelto — balanceo de carga en un clúster

Sistema de la GAA de la Sesión 6, con la matriz de acoplamiento entre tres microservicios:

```
 4x₁ +  2x₂ + 1x₃ = 14      (tráfico global)
12x₁ + 10x₂ + 5x₃ = 46      (procesamiento backend)
-8x₁ +  8x₂ + 7x₃ = 26      (consultas a base de datos)
```

Factorización por Doolittle, con det(A) = 48:

```
      1   0   0            4   2   1
L =   3   1   0      U =   0   4   2
     -2   3   1            0   0   3
```

| Vector de tráfico | b | y (progresiva) | x (solución) |
|---|---|---|---|
| b₁ — SLA inicial | (14, 46, 26) | (14, 4, 42) | (3, −6, 14) |
| b₂ — pico de tráfico | (20, 62, 30) | (20, 2, 64) | (4.75, −10.1667, 21.3333) |

El segundo vector no repite la descomposición: 54 flops con LU frente a 72 recalculando.
El x₂ negativo indica que el modelo pide liberar carga del microservicio S₂, señal de que
el SLA es inconsistente con los coeficientes de acoplamiento de ese nodo.
