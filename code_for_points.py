from pylab import *
import time

#####################################
###   PARAMETROS INICIALES   ########
#####################################

# Retícula de puntos.
Datos: list = []

# Tolerancia
tol: float = 1e-09
# Distancia radial
d: float = 0.01

#####################################
####     MÉTODO DE EULER   ##########
#####################################

# Número de Pasos
n: float = 2000
# Valor de trampa de escape
escape_radius: float = 100

#### LOCUS DE CONEXIDAD ####

# Número de iteraciones (órbita crítica)
num_iter: int = 200
# radio de escape del crítico
critic_radius: float = 100

# Punto de Inicial tomado del archivo soluciones_z0.py
pt_inicial: complex = -0.243965359233484 + 1.322407235713008 * 1j #z1_1
pt_initial1: complex = pt_inicial#-0.17708744616823316+1.3095100102028243*1j
pt_initial2: complex = pt_inicial#-0.38392300349284747+1.3354543400377912*1j
z0: complex = pt_initial1
w0: complex = pt_initial2

def cnan(z: complex):
    return isnan(z.real) or isnan(z.imag)

def F(x: complex, y: complex):
    f: complex = x*y
    g: complex = f + 1.0
    h: complex = x*x
    c: complex = f*g*(h+1.0)
    k: complex = c*g + 2.0
    return f*((c*k*g + 1.0)*(h*g*(h+1.0)) + k*c) + 1.0

def gradx(x:complex, y: complex):
    f: complex = x * y
    g: complex = f + 1.0
    h: complex = x * x
    j: complex = h + 1.0
    c: complex = f * g * j
    d: complex = h * g * j
    k: complex = c * g + 2.0
    e: complex = c * k * g
    c_x: complex = y * j * (g + f) + 2.0 * f * g * x
    d_x: complex = x * (2.0 * g * j + f * j + 2.0 * h * g)
    e_x: complex = c_x * k * g + c * g * (c_x * g + c * y) + c * k * y
    return e_x * d * f + e * d_x * f + e * d * y + c * f * (c_x * g + c * y) \
            + k * c_x * f + k * c * y

def grady(x: complex, y: complex):
    f: complex = x * y
    g: complex = f + 1.0
    h: complex = x * x
    j: complex = h + 1.0
    c: complex = f * g * j
    d: complex = h * g * j
    k: complex = c * g + 2.0
    e: complex = c * k * g
    c_y: complex = x * (g + f) * j
    d_y: complex = x * h * j
    e_y: complex = c_y * k * g + c * g * (c_y * g + c * x) + c * k * x
    return e_y * d * f + (e + 1.0) * d_y * f + (e + 1.0) * d * x \
        + (c_y * g + c * x) * c * f + k * c_y * f + k * c * x


## Funciones vector tangente
phi: complex = 1.0

def X1(x: complex, y: complex):
    return phi*gradx(x, y)

def X2(x: complex, y: complex):
    return -phi*grady(x, y)

### Función de Corrección (Newton-Raphson)
def correccion(z0: complex, w0: complex):
    z: complex = 0.0
    a: complex = gradx(z0, w0)
    b: complex = grady(z0, w0)
    i: int = 0
    while i < num_iter and abs(F(z0 + z*a, w0 + z*b)) > tol:
        z -= F(z0 + z*a, w0 + z*b)/(a*gradx(z0 + z*a, w0 + z*b) + b*grady(z0 + z*a, w0 + z*b))
        i += 1
    return z0 + z*a, w0 + z*b

# Método de Euler (integración) con correccion
def Euler(z0: complex, w0: complex, t0: float, t1: float, pasos: int, f1, f2):
    datos = []
    p0 = [z0, w0]
    datos.append(p0)
    h = (t1-t0)/pasos

    for i in range(pasos):
        a: complex = f1(p0[0],p0[1])
        b: complex = f2(p0[0], p0[1])
        p0 = [p0[0] + h*a, p0[1] + h*b]

        if cnan(p0[0]) or cnan(p0[1]):
            break
        if abs(F(p0[0], p0[1])) > escape_radius:
            break
        if abs(F(p0[0], p0[1])) > tol:
            p0 = correccion(p0[0], p0[1])
        datos.append(p0)
    return datos

def Euler_punto(z0: complex, w0: complex, t0: float, t1: float, pasos: int, f1, f2, theta, rad):
    datos = []
    z0 = z0 * exp(2.0*pi*1j*theta/n)
    p0 = [z0, w0]
    datos.append(p0)
    h = (t1-t0)/pasos

    for i in range(rad):
        a: complex = f1(p0[0],p0[1])
        b: complex = f2(p0[0], p0[1])
        p0 = [p0[0] + h*a, p0[1] + h*b]

        if cnan(p0[0]) or cnan(p0[1]):
            break
        if abs(F(p0[0], p0[1])) > escape_radius:
            break
        if abs(F(p0[0], p0[1])) > tol:
            p0 = correccion(p0[0], p0[1])
        datos.append(p0)
    return datos

# funcion de velocidades de escape
def velocidad_escape(x: complex, y: complex):
    if cnan(x) or cnan(y):
        return -1
    # cambio de coordenadas
    a: complex = -x/3.0 + y/3.0
    v: complex = 2.0*x/3.0 + y/3.0

    i: int = 0
    z: complex = -a
    while i < num_iter and abs(z) < critic_radius:
        i += 1
        z = ((z-a)**2)*(z + 2.0*a) + v
    return i

# main
# Programa principal
if __name__ == "__main__":
    inicio = time.time()
    file = open("txt_generados/test_new0.txt", "w")
    file1 = open("txt_generados/slice_points0.txt","w")
    # aumento de ángulo
    delta_phi: complex = exp((2.0/n)*pi*1j)
    # recorrer en dirección del nuevo phi y guardar velocidades de escape de cero  para los pols
    for i in range(n):
        datos = Euler(z0, w0, 0, d, n, X1, X2)
        j = 0
        N = len(datos)
        print("vamos en el paso ", i, " de ", n)
        while j < N-1:
            v = str(velocidad_escape(datos[j][0],datos[j][1]))
            file.write(v + " ")
            file1.write(str(datos[j][0]) + "," + str(datos[j][1]) + " ")
            j+=1
        while j < n:
            file.write("-1" + " ")
            file1.write("NaN" + " ")
            j+=1
        file.write("\n")
        file1.write("\n")
        phi *= delta_phi
        Datos.append(datos)
    file.close()
    fin = time.time()
    print(f"Tiempo de proceso: {(fin-inicio)/60:2f} minutos")