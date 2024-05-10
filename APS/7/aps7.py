import matplotlib.pyplot as plt
from typing import Any
from math import sqrt
from numpy import ndarray, floating, array, zeros, ix_, delete, dot, copy, inf, all
from numpy.linalg import solve, norm
import yaml


class No:

    def __init__(self, x: float, y: float) -> None:
        """
        Parâmetros:
            x (float): Coordenada x do nó
            y (float): Coordenada y do nó
        """
        self.x = x
        self.y = y


class Elemento:

    def __init__(self, comeco: No, fim: No, graus_liberdade: list[int], E: float, A: float, max_tracao: float, max_compressao: float) -> None:
        """
        Parâmetros:
            comeco (No): Nó de início do elemento
            fim (No): Nó de fim do elemento
            graus_liberdade (list[int]): Lista de graus de liberdade do elemento
            E (float): Módulo de elasticidade do material do elemento
            A (float): Área da seção transversal do elemento
            max_tracao (float): Tensão de ruptura a tração do elemento
            max_compressao (float): Tensão de ruptura a compressão do elemento
        """

        self.epsilon = 0
        self.sigma = 0

        self.nos = [comeco, fim]

        self.E = E
        self.A = A
        self.max_tracao = max_tracao
        self.max_compressao = max_compressao

        length = sqrt((fim.x - comeco.x)**2 + (fim.y - comeco.y)**2)
        cos = (fim.x - comeco.x)/length
        sen = (fim.y - comeco.y)/length

        self.length = length
        self.cos = cos
        self.sen = sen

        self.rigidez = E*A/length * array([[cos**2, cos*sen, -cos**2, -cos*sen],
                                           [cos*sen, sen**2, -cos*sen, -sen**2],
                                           [-cos**2, -cos*sen, cos**2, cos*sen],
                                           [-cos*sen, -sen**2, cos*sen, sen**2]])

        self.graus_liberdade = array(graus_liberdade)-1

    def verifica(self, U: ndarray[floating[Any]]) -> bool:
        """
        Método que verifica se o elemento está dentro dos limites de tensão

        Argumentos:
            U (ndarray[float]): Vetor de deslocamentos

        Retorna:
            bool: True se o elemento está dentro dos limites de tensão, False caso contrário
        """
        self.epsilon = 1/self.length * \
            dot(array([-self.cos, -self.sen, self.cos, self.sen]),
                U[ix_(self.graus_liberdade)])

        self.sigma = self.E*self.epsilon

        if self.sigma > self.max_tracao or self.sigma < -self.max_compressao:
            return False

        return True


class Trelica:

    def __init__(self, elementos: list[Elemento]) -> None:
        self.elementos = elementos

        self.K = zeros((graus_liberdade, graus_liberdade))

        for elemento in elementos:
            self.K[ix_(elemento.graus_liberdade,
                       elemento.graus_liberdade)] += elemento.rigidez

def gauss_seidel(K_restrito, P_restrito, max_de_iteracoes=1000, tolerancia=1e-10):
    U_restrito = zeros(len(P_restrito))
    
    for _ in range(max_de_iteracoes):
        U_novo = copy(U_restrito)
        
        for i in range(len(P_restrito)):
            somatoria_1 = dot(K_restrito[i, :i], U_novo[:i])
            somatoria_2 = dot(K_restrito[i, i+1:], U_restrito[i+1:])
            
            U_novo[i] = (P_restrito[i] - somatoria_1 - somatoria_2) / K_restrito[i, i]
        
        # Verifica se já convergiu
        if norm(U_novo - U_restrito, ord=inf) < tolerancia:
            return U_novo
        
        U_restrito = U_novo
    
    return U_restrito


nos: list[No] = []
elementos: list[Elemento] = []

with open("dados.yml", "r") as stream:
    data = yaml.safe_load(stream)

i = 1

for no in data["Nós"]:
    nos.append(No(no[0], no[1]))
    print(f"Adicionando nó {i} ({no[0]}, {no[1]})")
    i += 1

graus_liberdade = len(nos)*2

A = data["Área da Seção Transversal"]
E = data["Módulo de Elasticidade"]
MAX = data["Tensão Ruptura"]

print()
i = 1

for elemento in data["Elementos"]:
    comeco = elemento[0]
    fim = elemento[1]

    graus_de_liberdade_do_elemento = [comeco*2-1, comeco*2, fim*2-1, fim*2]

    elementos.append(
        Elemento(nos[comeco-1], nos[fim-1], graus_de_liberdade_do_elemento, E, A, MAX, MAX))

    print(f"Adicionando elemento {i} ({comeco}, {fim})")

    i += 1

P = zeros(graus_liberdade)

for forca in data["Forças"]:
    P[forca[0]] += forca[1]

restricoes = array(data["Restrições"])-1

K = zeros((graus_liberdade, graus_liberdade))

for elemento in elementos:
    K[ix_(elemento.graus_liberdade, elemento.graus_liberdade)] += elemento.rigidez

K_reduzido = delete(K, restricoes, axis=0)
K_reduzido = delete(K_reduzido, restricoes, axis=1)

P_reduzido = delete(P, restricoes)


graus_de_liberdade_sem_restricoes = [
    x for x in array(range(graus_liberdade))-1 if x not in restricoes]

graus_de_liberdade_sem_restricoes.remove(-1)

U = zeros(graus_liberdade)

#U[ix_(graus_de_liberdade_sem_restricoes)] += solve(K_reduzido, P_reduzido)
U[ix_(graus_de_liberdade_sem_restricoes)] += gauss_seidel(K_reduzido, P_reduzido)

print()

for s, g in zip(solve(K, P), gauss_seidel(K, P)):
    print(s, g)

P = dot(K, U)

print()
print("Deslocamentos:")

for i in range(len(nos)):
    print(f"O Deslocamento do nó {i+1} é: ({U[2*i-1]}, {U[2*i]})")

print()
print("Tensoes:")

for i in range(len(elementos)):
    fail = not elementos[i].verifica(U)

    sigma = elementos[i].sigma

    if sigma < 0:
        if fail:
            print(f"O elemento {
                  i+1} não suporta a compressão de {-sigma/1e6:.2f} MPa")
        else:
            print(f"O elemento {
                  i+1} está em compressão de {-sigma/1e6:.2f} MPa")
    else:
        if fail:
            print(f"O elemento {
                  i+1} não suporta a tração de {sigma/1e6:.2f} MPa")
        else:
            print(f"O elemento {i+1} está em tração de {sigma/1e6:.2f} MPa")

print()
print("Reações de apoio:")

for restricao in restricoes:
    print(f"A reação de apoio no grau de liberdade {
          restricao+1} é: {P[restricao]}")
