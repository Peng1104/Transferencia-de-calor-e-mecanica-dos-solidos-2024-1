import numpy as np

# Constates do problema
k = 230  # Condutividade térmica em W/(m·°C) (condução)
c = 897  # Calor específico em J/(kg·°C) (capacidade térmica)
rho = 2700  # Densidade em kg/m³
alpha = k / (c * rho)  # Difusividade térmica em m²/s

tempo = 100  # Tempo de simulação em segundos
dt = 1e-3 # Delta T

altura = 0.4
dy = altura/9 # Delta Y
ny = int(altura / dy) + 1

largura = 0.4
dx = largura/9 # Delta X
nx = int(largura / dx) + 1

# Matriz de temperatura
temperatura = np.zeros((ny, nx))  # Create a grid of nx by ny points set to 0°C

# Condições de contorno
temperatura[0, :] = 150  # Temperatura na borda superior é de 150°C
y_start = 1

temperatura[-1, :] = 0   # Temperatura na borda inferior é de 0°C
y_end = ny - 1

# temperatura[:, 0] = X  # Temperatura na borda esquerda é isolada
x_start = 0

temperatura[:, -1] = 50  # Temperatura na borda direita é de 50°C
x_end = nx - 1

# Simulação
for t in range(int(tempo / dt)):
    prox_temperatura = temperatura.copy()
    
    for i in range(y_start, y_end):
        for j in range(x_start, x_end):
            if j == 0:
                # Equação para borda isolada
                prox_temperatura[i, j] = temperatura[i, j] + alpha * dt / dx**2 * (
                    2 * temperatura[i, j] + temperatura[i, j + 1]
                    + temperatura[i, j - 1] - 4 * temperatura[i, j]
                )
            else:
                # Equação de difusão
                prox_temperatura[i, j] = temperatura[i, j] + alpha * dt / dx**2 * (
                    temperatura[i + 1, j] + temperatura[i - 1, j] + 
                    temperatura[i, j + 1] + temperatura[i, j - 1] - 4 * temperatura[i, j]
                )
    temperatura = prox_temperatura

np.set_printoptions(precision=8)
print(temperatura)