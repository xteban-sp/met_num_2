"""
Métodos Numéricos · Unidad 2 · Sesión 7
Método de Jacobi — balance de carga en un clúster de 4 servidores.

x_i : peticiones asignadas al servidor i (miles por minuto)
Vector inicial x(0) = [0, 0, 0, 0]
Criterio de parada: error relativo ||x(k) - x(k-1)||inf / ||x(k)||inf < 1e-4

Fernando Istaña Flores · Universidad Peruana Unión
"""


def es_diagonal_dominante(A):
    """Verifica |a_ii| > sum_{j != i} |a_ij| en cada fila y devuelve el detalle."""
    detalle = []
    ok = True
    for i, fila in enumerate(A):
        diag = abs(fila[i])
        resto = sum(abs(v) for j, v in enumerate(fila) if j != i)
        cumple = diag > resto
        ok = ok and cumple
        detalle.append((i + 1, diag, resto, cumple))
    return ok, detalle


def norma_jacobi(A):
    """||T_J||inf = max_i sum_{j != i} |a_ij| / |a_ii|  (cota de convergencia)."""
    return max(sum(abs(v) for j, v in enumerate(f) if j != i) / abs(f[i])
               for i, f in enumerate(A))


def jacobi(A, b, x0, tol=1e-4, max_iter=100):
    """
    Método de Jacobi. Cada componente nueva usa SOLO valores de la iteración anterior:
        x_i(k) = ( b_i - sum_{j != i} a_ij * x_j(k-1) ) / a_ii
    Devuelve la solución y el historial [(k, x(k), error)].
    """
    n = len(b)
    x_ant = list(x0)
    historial = [(0, list(x0), None)]
    for k in range(1, max_iter + 1):
        x_nuevo = [0.0] * n
        for i in range(n):
            suma = sum(A[i][j] * x_ant[j] for j in range(n) if j != i)
            x_nuevo[i] = (b[i] - suma) / A[i][i]
        dif = max(abs(x_nuevo[i] - x_ant[i]) for i in range(n))
        norma = max(abs(v) for v in x_nuevo)
        error = dif / norma if norma != 0 else dif
        historial.append((k, x_nuevo, error))
        if error < tol:
            return x_nuevo, historial
        x_ant = x_nuevo
    return x_ant, historial


if __name__ == "__main__":
    A = [[10, -2, -1,  0],
         [-1,  8,  0, -2],
         [-2,  0, 12, -3],
         [ 0, -1, -2,  9]]
    b = [15, 18, 25, 20]
    x0 = [0, 0, 0, 0]

    print("=" * 78)
    print("ANÁLISIS PRELIMINAR — ¿A es estrictamente diagonal dominante por filas?")
    print("=" * 78)
    ok, det = es_diagonal_dominante(A)
    for fila, diag, resto, cumple in det:
        print("  Fila %d: |a_ii| = %2d  >  suma resto = %d   %s"
              % (fila, diag, resto, "cumple" if cumple else "NO cumple"))
    print("  Resultado: A %s EDD" % ("ES" if ok else "NO es"))
    print("  ||T_J||inf = %.4f < 1  ->  Jacobi converge para cualquier x(0)" % norma_jacobi(A))

    print()
    print("=" * 78)
    print("ITERACIONES DE JACOBI  (tol = 1e-4, error relativo en norma infinito)")
    print("=" * 78)
    x, hist = jacobi(A, b, x0, tol=1e-4)
    print("  k        x1          x2          x3          x4          error")
    for k, xk, err in hist:
        e = "     —" if err is None else "%.6e" % err
        print("%3d  %10.6f  %10.6f  %10.6f  %10.6f   %s" % (k, xk[0], xk[1], xk[2], xk[3], e))

    print()
    print("Solución aproximada tras %d iteraciones:" % (len(hist) - 1))
    for i, v in enumerate(x):
        print("  x%d = %.6f  miles de peticiones/min" % (i + 1, v))
    r = [sum(A[i][j] * x[j] for j in range(4)) - b[i] for i in range(4)]
    print("  residuo A·x − b = (" + ", ".join("%.2e" % v for v in r) + ")")
