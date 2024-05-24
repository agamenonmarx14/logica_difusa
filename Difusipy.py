import math
import scipy.integrate as integrate

class Fuzzifier:
    def cuadratica(l, x):
        b=l[0]
        a=l[1]
        if (b-a)<=x<=(b+a):
            u=(1-pow(((x-b)/a), 2))
        else:
            u=0
        return u

    def cuadraticaL(l, x):
        b=l[0]
        a=l[1]
        if x<b:
            u=1
        elif b<=x<=(b+a):
            u=(1-pow(((x-b)/a), 2))
        elif x>(b+a):
            u=0
        return u

    def cuadraticaR(l, x):
        b=l[0]
        a=l[1]
        if x<(b-a):
            u=0
        elif (b-a)<=x<=b:
            u = u=(1-pow(((x-b)/a), 2))
        elif x>b:
            u = 1
        return u


    def gaussiana(l, x):
        """b: valor para el que y=1
        xk: punto en x para el que y:0.5"""
        a=abs(l[1]-l[0])/(math.sqrt(math.log(2)))
        u=math.exp(-((x-l[0])/a)**2)
        return u

    def sigmoide(l, x):
        """b=valor para el que y=0.99
        xk: punto en x para el que y:0.5"""
        a=math.log(99)/(l[0]-l[1])
        try:
            d=(1+math.exp(-a*(x-l[1])))
            u=1/d
        except:
            u=1
        return u


    def trapecial(l, x):
        """x: Valor de variable a evaluar
           a: Punto a la izquierda en x donde y=0
           b: Punto a la izquierda en x donde y=1
           c: Punto a la derecha en x donde y=1
           d: Punto a la derecha en x donde y=0
        """
        a=l[0]
        b=l[1]
        c=l[2]
        d=l[3]
        if (a<=x<b):
            u = (x-a)/(b-a)
        elif (b<=x<c):
            u = 1
        elif (c<=x<=d):
            u = (d-x)/(d-c)
        else:
            u = 0
        return u

    def triangular(l, x):
        """x: Valor de la variable a evaluar
           a: Ancho de la bisectriz en la base
           e: Punto en x de la bisectriz
           Retorna: Membresía en el rango dado
           """
        a=l[0]
        e=l[1]
        if (e-a)<=x<(e+a):
            u = (a-abs(x-e))/a
        else:
            u = 0
        return u

    def extremo(l, x):
        """a0: Punto en x para el que y=0
           a1: Punto en x para el que y=1
           Retorna: Membresía en el rango dado"""
        a0=l[0]
        a1=l[1]
        if (a0<a1):
            if (a1<=x):
                u=1
            elif (a0<x<a1):
                u=(x-a0)/(a1-a0)
            else:
                u=0
        if (a1<a0):
            if (x<=a1):
                u=1
            elif (a1<x<a0):
                u=(a0-x)/(a0-a1)
            else:
                u=0
        return u


class Operador_TNorma:
    def operadorMin(l):
        """l: Lista de valores a evaluar"""
        minvalor = min(l)
        return minvalor

    def prodHamacher(l):
        #Recibe una lista de funciones de membresía y devuelve T-Norma correspondiente.
        r=0
        for i in range (len(l)):
            if (i<1):
                r=(l[i]*l[i+1])/(l[i]+l[i+1]-l[i]*l[i+1])
            elif 2<=i:
                r=(r*l[i])/(r+l[i]-r*l[i])
            else:
                pass
        return r

    def prodEinstein(l):
        r=0
        for i in range (len(l)):
            if (i<1):
                r=(l[i]*l[i+1])/(2-(l[i]+l[i+1]-l[i]*l[i+1]))
            elif 2<=i:
                r=(r*l[i])/(2-(r+l[i]-r*l[i]))
            else:
                pass
        return r

    def producto(l):
        r=0
        for i in range (len(l)):
            if (i<1):
                r=l[i]*l[i+1]
            elif 2<=i:
                r=r*l[i]
            else:
                pass
        return r

class Operador_SNorma:
    def operadorMax(l):
        """l: Lista de valores a evaluar"""
        maxvalor = max(l)
        return maxvalor

    def suma_Hamacher(l):
        c=0
        for i in range(len(l)):
            if (i<1):
                c=(l[i]+l[i+1] - 2*l[i]*l[i+1])/(1-l[i]*l[i+1])
            elif 2<=i:
                c=(c+l[i] - 2*c*l[i])/(1-c*l[i])
            else:
                pass
        return c

    def suma_Einstein(l):
        c=0
        for i in range(len(l)):
            if (i<1):
                c=(l[i]+l[i+1])/(1+l[i]*l[i+1])
            elif 2<=i:
                c=(c+l[i])/(1+l[i]*c)
            else:
                pass
        return c

class Defuzzifier:
    def defgaussiano(l, u):
        b=l[0]
        xk=l[1]
        a=abs(l[1]-l[0])/(math.sqrt(math.log(2)))
        Ax=integrate.quad(lambda x: x*math.exp(-((x-b)/a)**2), xk, (2*b-xk))
        A=integrate.quad(lambda x: math.exp(-((x-b)/a)**2), xk, (2*b-xk))
        return [Ax[0]*u, abs(A[0]*u)]

    def defsigmoide(l, u):
        b=l[0]
        xk=l[1]
        a=math.log(99)/(b-xk)
        if xk<b:
            Ax=integrate.quad(lambda x: x/(1+math.exp(-a*(x-b))), xk, b)
            A=integrate.quad(lambda x: 1/(1+math.exp(-a*(x-b))), xk, b)
        else:
            Ax=integrate.quad(lambda x: x/(1+math.exp(-a*(x-b))), b, xk)
            A=integrate.quad(lambda x: 1/(1+math.exp(-a*(x-b))), b, xk)
        return [Ax[0]*u, abs(A[0]*u)]

    def cuadratica(l, u):
        b=l[0]
        xk=l[1]
        if xk<b:
            a=(b-xk)*math.sqrt(2)
            Ax=integrate.quad(lambda x: x*(1-pow(((x-b)/a), 2)), xk, (2*b-xk))
            A=integrate.quad(lambda x: (1-pow(((x-b)/a), 2)), xk, (2*b-xk))
        else:
            a=(xk-b)*math.sqrt(2)
            Ax=integrate.quad(lambda x: x*(1-pow(((x-b)/a), 2)), (2*b-xk), xk)
            A=integrate.quad(lambda x: (1-pow(((x-b)/a), 2)), (2*xk-b), xk)
        return [Ax[0]*u, A[0]*u]

    def cuadraticaL(l, u):
        b=l[0]
        xk=l[1]
        a=(xk-b)*math.sqrt(2)
        Ax=integrate.quad(lambda x: x*(1-pow(((x-b)/a), 2)), b, xk)
        A=integrate.quad(lambda x: (1-pow(((x-b)/a), 2)), b, xk)
        return [Ax[0]*u, A[0]*u]

    def cuadraticaR(l, u):
        b=l[0]
        xk=l[1]
        a=(b-xk)*math.sqrt(2)
        Ax=integrate.quad(lambda x: x*(1-pow(((x-b)/a), 2)), xk, b)
        A=integrate.quad(lambda x: (1-pow(((x-b)/a), 2)), xk, b)
        return [Ax[0]*u, A[0]*u]

    def triangular(l, u):
        b=l[0]
        xk=l[1]
        A=0.25*(5*b-3*xk)*u
        Ax=0.25*(5*b*b-3*b*xk)*u
        return [Ax, A]

    def semitriangular(l, u):
        b=l[0]
        xk=l[1]
        if b<xk:
            A = 0.75*(xk-b)*u
            Ax = 0.75*(b*xk-b*b)*u
        else:
            A =  0.75*(b-xk)*u
            Ax = 0.75*(b*b - b*xk)*u
        return [Ax, A]

class Operador_Generico:
    def operadorCON(l):
        """l: Lista de valores a evaluar"""
        con = []
        for i in l:
            con.append(i*i)
        return con

    def operadorDIL(l):
        """l: Lista de valores a evaluar"""
        dil = []
        for i in l:
            dil.append(pow(i, 0.5))
        return dil

class Regla:
    def regla(u):
        if (0.01<=u):
            j=True
        else:
            j=False
        return j
