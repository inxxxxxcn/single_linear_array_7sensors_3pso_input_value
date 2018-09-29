import matplotlib.pyplot as plt
import numpy as np
import sympy as sym


def b(x0, y0, z0, x1, y1, z1, m, I, D):
    # 目标到传感器的距离
    x, y, z = x1 - x0, y1 - y0, z1 - z0
    # 系数
    u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
    p = u0 / (4 * np.pi)
    # 磁矩大小和方向
    m = 152
    j, k, l = sym.cos(I) * sym.cos(D), sym.cos(I) * sym.sin(D), sym.sin(I)
    # 地磁方向
    I0, D0 = 1.1, -0.18
    u, v, w = sym.cos(I0) * sym.cos(D0), sym.cos(I0) * sym.sin(D0), sym.sin(I0)
    # 参数变形
    a = 2 * j * u - k * v - l * w
    b = 2 * k * v - j * u - l * w
    c = 2 * l * w - j * u - k * v
    d = 3 * (k * u + j * v)
    e = 3 * (l * u + j * w)
    f = 3 * (l * v + k * w)
    g = j * u + k * v + l * w
    r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
    # ba = p * m * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * (
    #         x ** 2 + y ** 2 + z ** 2) ** (-5 / 2)
    # ba = p * m * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-5)
    # ba = p * m * r ** (-3) * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-2)
    # 范围(1,2)
    btr = (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z)
    # 磁异常
    ba = p * m * r ** (-5) * btr
    return ba


if __name__ == "__main__":
    x1, y1, z1 = sym.symbols("x1, y1, z1")
    x0, m, I, D = sym.symbols("x0, m, I, D")
    ba = b(x0, 0, 0, x1, y1, z1, m, I, D)[0]
    ba = ba.subs(I, 2.13).subs(D, 1.28)
    ba = ba.subs(x1, 200).subs(y1, 200).subs(z1, 0)
    print(ba)
    # 画图
    size = 10
    xi = np.linspace(-50, 50, size)
    Gp = np.zeros(shape=size-1)
    for i in range(size - 1):
        Gp[i] = ba.subs(x0, xi[i]) - ba.subs(x0, xi[i + 1])
        print(Gp[i])
    xi = np.delete(xi, size-1, axis=0)
    plt.plot(xi, Gp)
    plt.show()
