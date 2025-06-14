#include <iostream>
#include <vector>
#include <complex>
#include <cmath>    // Operaciones
#include <iomanip>  // Configurar el número de decimales con el que se imprimen los resultados
#include <fstream>  // Abrir y cerrar archivos
using namespace std;
using complexd = complex<double>;
using tupla = pair<complexd, complexd>;

// Unidad imaginaria:
const complexd I (0,1);
// Pi:
const double pi = acos(-1);
// Tolerancia:
const double tol = 1e-09;
// Distancia:
const double d = 0.01;

// Método de Euler:

// Número de pasos de Euler:
const int n = 2000;
// Si el valor absoluto de un punto es mayor a este valor (i.e. se alejó mucho de la curva)
// detenemos el loop for de Euler:
const double demasiado_grande = 5;

// Determinar si un punto pertenece al locus de conexidad:

// Número de iteraciones:
const int num_iter = 200;
// Si la imagen de -a bajo un iterado del polinomio es mayor a este valor, 
// asumimos que -a tiende a infinito:
const double A = 100;

// Punto de partida: 
const complexd z0 = -0.243965359233484 + 1.322407235713008 * I;
const complexd w0 = z0;                                                                                  

bool cnan (complexd z) {
    return isnan(z.real()) || isnan(z.imag());
}

complexd F (complexd x, complexd y) {
    complexd f = x*y;
    complexd g = f + 1.0;
    complexd h = x*x;
    complexd c = f * g * (h + 1.0);
    complexd k = c * g + 2.0; 
    return f * ((c * k * g + 1.0) * (h * g * (h + 1.0)) + k * c) + 1.0;
}

// Funciones coordenadas del gradiente de F(x,y):

complexd grad1 (complexd x, complexd y) {
    complexd f = x*y;
    complexd g = f + 1.0;
    complexd h = x*x;
    complexd j = h + 1.0;
    complexd c = f * g * j;
    complexd d = h * g * j;
    complexd k = c * g + 2.0; 
    complexd e = c * k * g;
    complexd c_x = y * j * (g + f) + 2.0 * f * g * x;
    complexd d_x = x * (2.0 * g * j + f * j + 2.0 * h * g);
    complexd e_x = c_x * k * g + c * g * (c_x * g + c * y) + c * k * y;
    return e_x * d * f + e * d_x * f + e * d * y + c * f * (c_x * g + c * y)
                                                    + k * c_x * f + k * c * y;
}

complexd grad2 (complexd x, complexd y) {
    complexd f = x*y;
    complexd g = f + 1.0;
    complexd h = x*x;
    complexd j = h + 1.0;
    complexd c = f * g * j;
    complexd d = h * g * j;
    complexd k = c * g + 2.0; 
    complexd e = c * k * g;
    complexd c_y = x * (g + f) * j;
    complexd d_y = x * h * j;
    complexd e_y = c_y * k * g + c * g * (c_y * g + c * x) + c * k * x;
    return e_y * d * f + (e + 1.0) * d_y * f + (e + 1.0) * d * x
            + (c_y * g + c * x) * c * f + k * c_y * f + k * c * x;
}

// Funciones coordenadas del vector tangente:

complexd phi = 1.0; // Ángulo en que rotamos el capo original

complexd X1 (complexd a, complexd b) {
    return phi * grad2(a, b);
}

complexd X2 (complexd a, complexd b) {
    return -phi * grad1(a, b);
}

// Corrección con Newton-Rapson:

tupla correccion (complexd z0, complexd w0) {
    complexd z = 0;
    complexd a = grad1(z0, w0);
    complexd b = grad2(z0, w0);
    int i = 0;
    while (abs(F(z0+z*a, w0+z*b)) > tol)
    {
        z -= F(z0+z*a, w0+z*b)
                / (a*grad1(z0+z*a, w0+z*b) + b*grad2(z0+z*a, w0+z*b)); 
        ++i;
    }
    return tupla (z0+z*a, w0+z*b);
}

// Implementación del método de Euler (con corrección):

vector<tupla> Euler (complexd z0, complexd w0, double t0, double t1, int pasos, 
                    complexd f1 (complexd a, complexd b), complexd f2 (complexd a, complexd b)) {

    vector<tupla > datos;
    tupla p0 (z0, w0);
    datos.push_back(p0);

    double h = (t1-t0)/pasos;

    for (int i = 1; i<=pasos; ++i)
    {
        complexd a = f1(p0.first, p0.second);
        complexd b = f2(p0.first, p0.second);
        p0.first += h*a;    // <-- Puede ser h*a/|a| : esto equivale a usar el campo tangente normalizado
        p0.second += h*b;

        // Corrección:
        if (cnan(p0.first) || cnan(p0.second))
        {
            // cout << "nan!" << endl;
            break;
        }
        if (abs(F(p0.first, p0.second)) > demasiado_grande)
            break;
        if (abs(F(p0.first, p0.second)) > tol)
            p0 = correccion(p0.first, p0.second);
        
        datos.push_back(p0);
        // cout << i << " " << abs(F(p0.first, p0.second)) << endl;
    }

    return datos;
}

// Función para determinar qué puntos pertenecen al locus de conexidad.

int velocidad_escape (complexd x, complexd y) {

    // Si el punto tiene coordenadas nan (i.e. se fue a infinito): 
    if (cnan(x) || cnan(y))
        return -1;

    // Cambio de coordenadas:
    complexd a = -x/3.0 + y/3.0;
    complexd v = 2.0*x/3.0 + y/3.0;

    // Aplicación iterada del polinomio:
    int i = 0;
    complexd z = -a;
    while (i<num_iter && abs(z) < A)
    {
        ++i;
        z = pow(z-a,2)*(z+2.0*a) + v;
    }

    return i;

}

// Programa principal

int main() {

    // Abrir archivo:
    ofstream myfile;
    myfile.open ("S4.txt");

    // Aumento ángulo:
    complexd delta_phi = exp((2.0/n) * pi * I);

    // Recorrer en la otra y guardar velocidades de escape:
    for (int i=0; i<=n; i += 1)
    {
        vector<tupla> datos = Euler (z0, w0, 0, d, n, X1, X2);
        // cout << "Columna " << i << endl;

        int j = 0;
        while (j < n + 1)
        {
            myfile << velocidad_escape(datos[j].first, datos[j].second) << " ";
            ++j;
        }
        while (j < n + 1)
        {
            myfile << -1 << " ";
            ++j;
        }
        myfile << "\n";

        phi *= delta_phi;
    }

    // Cerrar:
    myfile.close();
    
}
