import math


def calculate_e(value: float, min_error: float) -> tuple[float, int]:
    error = 1.0
    e = 1.0

    counter = 1

    while error > min_error:
        new_e = e + value**counter/math.factorial(counter)
        error = abs((e - new_e) / new_e)
        e = new_e
        counter += 1

    return e, counter - 1


if __name__ == "__main__":
    e, counter = calculate_e(1, 0.0005)

    print("The math library value of e is: ", math.e)
    print("The value of e is: ", e)
    print("The number of iterations is: ", counter)
