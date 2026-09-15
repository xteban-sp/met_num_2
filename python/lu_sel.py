"""
Métodos Numéricos · Unidad 2 · Sesión 6
Solución de sistemas de ecuaciones lineales por factorización LU.

Caso: balanceo de carga en un clúster de 3 microservicios.
    A x = b     A: matriz de acoplamiento    b: SLA de latencia    x: tiempos de cómputo

Métodos incluidos: Doolittle, Crout y Cholesky.
Fernando Istaña Flores · Universidad Peruana Unión
"""

import math


# --------------------------------------------------------------------------
# Factorizaciones
# --------------------------------------------------------------------------
def doolittle(A):
    """A = L U con L triangular inferior de diagonal unitaria."""
    n = len(A)
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [row[:] for row in A]
    for k in range(n):
        if abs(U[k][k]) < 1e-12:
            raise ValueError("pivote nulo en U[%d][%d]: se requiere pivoteo parcial" % (k, k))
        for i in range(k + 1, n):
            m = U[i][k] / U[k][k]      # multiplicador -> entrada de L
            L[i][k] = m
            for j in range(k, n):
                U[i][j] -= m * U[k][j]
            U[i][k] = 0.0
    return L, U


def crout(A):
    """A = L U con U triangular superior de diagonal unitaria."""
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for j in range(n):
        for i in range(j, n):
            L[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(j))
        if abs(L[j][j]) < 1e-12:
            raise ValueError("L[%d][%d] = 0: Crout no puede continuar" % (j, j))
        for i in range(j + 1, n):
            U[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(j))) / L[j][j]
    return L, U


def cholesky(A):
    """A = L Lt. Requiere A simétrica y definida positiva."""
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for j in range(n):
        rad = A[j][j] - sum(L[j][k] ** 2 for k in range(j))
        if rad <= 0:
            raise ValueError("radicando %.4f <= 0: la matriz no es definida positiva" % rad)
        L[j][j] = math.sqrt(rad)
        for i in range(j + 1, n):
            L[i][j] = (A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))) / L[j][j]
    U = [[L[j][i] for j in range(n)] for i in range(n)]
    return L, U


# --------------------------------------------------------------------------
# Sustituciones — sirven para los tres métodos
# --------------------------------------------------------------------------
def resolver(L, U, b):
    n = len(b)
    y = [0.0] * n                                   # L y = b  (progresiva)
    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
    x = [0.0] * n                                   # U x = y  (regresiva)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return y, x


def residuo(A, x, b):
    n = len(b)
    return [sum(A[i][j] * x[j] for j in range(n)) - b[i] for i in range(n)]


def mostrar(nombre, M):
    print("%s =" % nombre)
    for fila in M:
        print("   [" + "  ".join("%9.4f" % v for v in fila) + "]")


def vec(v):
    return "(" + ", ".join("%.4f" % c for c in v) + ")"


# --------------------------------------------------------------------------
# Caso de la guía
# --------------------------------------------------------------------------
if __name__ == "__main__":
    A = [[4, 2, 1],
         [12, 10, 5],
         [-8, 8, 7]]
    b1 = [14, 46, 26]        # SLA inicial
    b2 = [20, 62, 30]        # pico de tráfico

    print("=" * 66)
    print("FASE 1 — Descomposición A = LU por Doolittle (una sola vez)")
    print("=" * 66)
    L, U = doolittle(A)
    mostrar("L", L)
    mostrar("U", U)
    det = 1.0
    for i in range(len(A)):
        det *= U[i][i]
    print("det(A) = %.4f" % det)

    print()
    print("=" * 66)
    print("FASE 1 y 2 — Sustituciones para cada vector de tráfico")
    print("=" * 66)
    for nombre, b in [("b1", b1), ("b2", b2)]:
        y, x = resolver(L, U, b)
        r = residuo(A, x, b)
        print("%s = %s" % (nombre, vec([float(v) for v in b])))
        print("   y = %s      (sustitución progresiva)" % vec(y))
        print("   x = %s      (sustitución regresiva)" % vec(x))
        print("   residuo A·x − b = %s" % vec(r))
        print()

    print("=" * 66)
    print("COMPROBACIÓN — otros métodos de factorización")
    print("=" * 66)
    Lc, Uc = crout(A)
    yc, xc = resolver(Lc, Uc, b1)
    print("Crout    x = %s" % vec(xc))

    S = [[4, -2, 2],
         [-2, 10, -7],
         [2, -7, 30]]
    bs = [14, -43, 106]
    Ls, Us = cholesky(S)
    ys, xs = resolver(Ls, Us, bs)
    mostrar("L de Cholesky (matriz simétrica de ejemplo)", Ls)
    print("Cholesky x = %s" % vec(xs))

    print()
    n = len(A)
    fact = round(2 * n ** 3 / 3)
    sust = 2 * n ** 2
    print("Costo con LU        : %d + 2×%d = %d flops" % (fact, sust, fact + 2 * sust))
    print("Costo recalculando  : 2×(%d + %d) = %d flops" % (fact, sust, 2 * (fact + sust)))
