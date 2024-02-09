import math


def calculate_sin(value: float, min_error: float) -> tuple[float, int]:
    error = 1.0
    sin = 0.0

    counter = 1

    while error > min_error:
        new_sin = sin + ((-1)**(counter - 1) * value **
                         (2*counter - 1))/math.factorial(2*counter - 1)

        error = abs((sin - new_sin) / new_sin)
        sin = new_sin
        counter += 1

    return sin, counter - 1


if __name__ == "__main__":
    sin_x, counter = calculate_sin(math.pi/3, 5e-6)

    print("The math library value of sin(pi/3) is: ", math.sin(math.pi/3))
    print("The value of sin(pi/3) is: ", sin_x)
    print("The number of iterations is: ", counter)
