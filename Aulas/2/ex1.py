class Pair:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y


def newton(target: float, values: list[Pair]) -> float:
    n = len(values)

    coefficients: list[float] = []

    for i in range(n):
        coefficients.append(values[i].y)

        for j in range(i):
            coefficients[i] = (coefficients[i] -
                               coefficients[j]) / (values[i].x - values[j].x)

    result = coefficients[0]

    for i in range(1, n):
        term = 1

        for j in range(i):
            term *= (target - values[j].x)

        result += coefficients[i] * term

    return result


def lagrange(target: float, values: list[Pair]) -> float:
    n = len(values)

    result = 0.0

    for i in range(n):
        f_x = values[i].y

        for j in range(n):
            if j != i:
                result += f_x * (target - values[j].x) / (values[i].x - values[j].x)
    
    return result


# Volume procurado
volume = 1.34567

# Dados
valores = [Pair(1.31623, 7.78941), Pair(1.54934, 8.2236)]

# Interpolação de Newton
print(
    f"Volume calculado pela interpolação de Newton: {newton(volume, valores)}")

# Interpolação de Lagrange
print(
    f"Volume calculado pela interpolação de Lagrange: {lagrange(volume, valores)}")
