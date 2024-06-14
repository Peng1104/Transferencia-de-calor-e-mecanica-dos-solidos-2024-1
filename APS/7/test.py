from typing import Any
from numpy import ndarray, floating, zeros, dot, copy, array
from numpy.linalg import norm


def gauss_seidel(K_restrito, P_restrito, max_de_iteracoes=10000, tolerancia=1e-10) -> ndarray[floating[Any]]:
    """
    Método de Gauss-Seidel para resolver sistemas lineares
    Parâmetros:
        K_restrito (ndarray): Matriz de rigidez restrita
        P_restrito (ndarray): Vetor de forças restrito
        max_de_iteracoes (int): Número máximo de iterações
        tolerancia (float): Tolerância para convergência
    Retorna:
        ndarray: Vetor de deslocamentos restrito
    """

    size = len(P_restrito)
    U_restrito = zeros(size)

    for _ in range(max_de_iteracoes):
        U_antigo = copy(U_restrito)

        for i in range(size):
            # Calcula a somatória
            # K[i, j] * U[j] para j != i
            somatoria = dot(K_restrito[i, :], U_antigo) - \
                K_restrito[i, i] * U_restrito[i]  # remove a diagonal (j == i)

            U_restrito[i] = (P_restrito[i] - somatoria) / K_restrito[i, i]

        norma = norm(U_restrito)

        # Verifica se já convergiu
        if norma > 0:
            if norm(U_restrito - U_antigo) / norma < tolerancia:
                return U_restrito

    return U_restrito


M = array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
V = array([7.85, -19.3, 71.4])

print(gauss_seidel(M, V))
