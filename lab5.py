import numpy as np
import matplotlib.pyplot as plt

def graph(a, b, c, d, arg, x0):
	return a + b * (arg - x0) + c * (arg - x0) ** 2 + d * (arg - x0) ** 3

def coefC(cArr, mass, n, h):
	A = [0] * (n - 2)
	for i in range(n - 2):
		A[i] = [0] * (n - 2)
	for i in range(1, n - 3):
		j = i - 1
		A[i][j] = h
		A[i][j + 1] = 4 * h
		A[i][j + 2] = h
	A[0][0] = 4 * h
	A[0][1] = h
	A[n - 3][n - 4] = h
	A[n - 3][n - 3] = 4 * h

	F = [0] * (n - 2)
	for i in range(1, n - 1):
		F[i - 1] = 3 * ((mass[i + 1][1] - mass[i][1]) / h - (mass[i][1] - mass[i - 1][1]) / h)
	for i in range(1, n - 2):
		F[i] -= (A[i][i - 1] * F[i - 1] / A[i - 1][i - 1])
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1])
	for i in range(1, n - 2):
		F[i] -= (A[i][i - 1] * F[i - 1] / A[i - 1][i - 1])
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1])
        
	cArr[n - 2] = F[n - 3] / A[n - 3][n - 3]
	for i in range(n - 4, -1, -1):
		cArr[i + 1] = (F[i] - A[i][i + 1] * cArr[i + 2]) / A[i][i]
	cArr[0] = 0
	cArr[n - 1] = 0

def splineFunc(arg, a, b, c, d, func):
	i = 0
	for j in range(1, len(c) - 1):
		if (arg > func[j][0]):
			i += 1
	return a[i] + b[i] * (arg - func[i][0]) + c[i] * (arg - func[i][0]) ** 2 + d[i] * (arg - func[i][0]) ** 3

def trapeze(x0, xn, step, func, a, b, c, d):
    sum = 0
    delta = (xn - x0) / step
    while(x0 <= xn):
        sum += (splineFunc(x0, a, b, c, d, func) + splineFunc(x0 + delta, a, b, c, d, func)) * delta / 2
        x0 += delta
    return sum

def simpson(x0, xn, step, func, a, b, c, d):
    sum = 0
    delta = (xn - x0) / step
    while (x0 <= xn):
        sum += delta * (splineFunc(x0, a, b, c, d, func) + 4 * splineFunc(x0 + delta / 2, a, b, c, d, func) + splineFunc(x0 + delta, a, b, c, d, func)) / 6
        x0 += delta
    return sum

def rectangle(x0, xn, step, func, a, b, c, d):
	sum = 0
	delta = (xn - x0) / step
	while (x0 <= xn):
		sum += splineFunc(x0 + delta / 2, a, b, c, d, func) * delta
		x0 += delta
	return sum

def runge(meth, x0, xn, step, func, a, b, c, d):
	if (meth == 1):
		return np.fabs(trapeze(x0, xn, step, func, a, b, c, d) - trapeze(x0, xn, 2 * step, func, a, b, c, d)) / 3
	elif (meth == 2):
		return np.fabs(simpson(x0, xn, step, func, a, b, c, d) - simpson(x0, xn, 2 * step, func, a, b, c, d)) / 15
	elif (meth == 3):
		return np.fabs(rectangle(x0, xn, step, func, a, b, c, d) - rectangle(x0, xn, 2 * step, func, a, b, c, d)) / 3

n = 9
h = 0.25
eps = 1e-4
func = [0] * n
for i in range (n):
    func[i] = [0] * 2

j = 0
for i in range(0, n):
    func[i][0] = j
    j += h
func[0][1] = 1
func[1][1] = 0.979915
func[2][1] = 0.927295
func[3][1] = 0.858001
func[4][1] = 0.785398
for i in range(5, n):
	func[i][1] = 0.716844

a = [0] * (n - 1)
b = [0] * (n - 1)
c = [0] * n
d = [0] * (n - 1)

coefC(c, func, n, h)
for i in  range(n - 1):
	a[i] = func[i][1]
	b[i] = (func[i + 1][1] - func[i][1]) / h - c[i] * h - (c[i + 1] - c[i]) * h / 3
	d[i] = (c[i + 1] - c[i]) / (3 * h)

print("Trapeze method:")
integ = 0
pogr = 1
step = 50
while (pogr > eps):
	integ = trapeze(func[0][0], func[n - 1][0], step, func, a, b, c, d)
	pogr = runge(1, func[0][0], func[n - 1][0], step, func, a, b, c, d)
	step *= 2
print(integ, pogr)

print("Simpson method:")
integ = 0
pogr = 1
step = 50
while (pogr > eps):
	integ = simpson(func[0][0], func[n - 1][0], step, func, a, b, c, d)
	pogr = runge(2, func[0][0], func[n - 1][0], step, func, a, b, c, d)
	step *= 2
print(integ, pogr)
 
print("Rectangle method:")
integ = 0
pogr = 1
step = 50
while (pogr > eps):
	integ = rectangle(func[0][0], func[n - 1][0], step, func, a, b, c, d)
	pogr = runge(3, func[0][0], func[n - 1][0], step, func, a, b, c, d)
	step *= 2
print(integ, pogr)

for i in range(n - 1):
    arg = np.linspace(func[i][0], func[i+1][0], 100)
    plt.plot(arg, graph(a[i], b[i], c[i], d[i], arg, func[i][0]), 'g')
for i in range(n):
    plt.plot(func[i][0],func[i][1],'r*')
plt.grid()
plt.show()