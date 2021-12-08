#Библиотеки:
from sys import stdin
#Импортируем numpy для поиска точного решения:
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Итерационная формула для квадратного корня:
def root(a, eps = 10 ** (-8), w0 = 0):
    if (a>=0):
        w1 = abs(eps) + abs(w0) + 10 ** (-8)
        while (abs(w1 - w0)>abs(eps)):
            w0 = w1
            w1 = 0.5*(w0 + a/w0)
    else:
        print("Negative value")
    return w1
#
                                                            #СОЗДАНИЕ КЛАССА МАТРИЦ:

class Matrix():
    
    #Конструктор класса
    def __init__(self1, A):
        self1.matrix = A.copy()
        self1.nrows = len(self1.matrix)
        self1.ncols = len(self1.matrix[0])
    #


    #Единичная матрица. Или матрица с d на диагонали.
    @staticmethod
    def eye(n, d = 1):

        E = [[0] * n for i in range(n)]
        for i in range(n):
            for j in range(n): 
                if (i==j):
                    E[i][j] = float(d)
        return Matrix(E)
    #



    #Перестановка строк
    def swap(self2, index1, index2):

        A = list(self2.matrix).copy()
        t = A[index1].copy()
        A[index1] = list(A[index2]).copy()
        A[index2] = list(t).copy()

        return Matrix(A)
    #


    #Скалярное умножение строк
    @staticmethod
    def scal_mul(row1, row2):
        s = 0

        if (len(row1) == len(row2)):
            for i in range(len(row1)):
                s = s + row1[i]*row2[i]
        else:
            print("Wrong dimensions!")
        return s
    #

    #Евклидова норма вектора:
    @staticmethod
    def norm_vect(row1):
        return root(Matrix.scal_mul(row1, row1))
    
    #

    #Транспонирование:
    def trans(self3):
        return [list(i) for i in zip(*self3.matrix)]
    #

    #Максимум в списке:
    def maxim(a):
        t = a[0]
        for i in range(len(a)):
            if (a[i] > t):
                t = a[i]
        return t
    #

    #Максимум в списке:
    def minim(a):
        t = a[0]
        for i in range(len(a)):
            if (a[i] < t):
                t = a[i]
        return t
    #

    #Сумма всех элементов:
    def summa(a):
        s = 0
        for i in range(len(a)):
                s = s + a[i]
        return s
    #

    #Модуль всех элементов:
    def mdl(A):
        return [abs(a) for a in A]
    #

    #Норма ||A||1 (максимальная сумма столбцов):
    def Norma_1(A):
        return Matrix.maxim([Matrix.summa(a) for a in [[abs(a[i]) for a in A.matrix] for i in range(A.ncols)]])
    #

    #
    #


    #Норма ||A||inf (максимальная сумма строк):
    def Norma_2(A):
        return Matrix.maxim([Matrix.summa(Matrix.mdl(a)) for a in A.matrix])
    #

    #Евклидова норма ||A||2:
    def Norma_3(A):
        s = 0
        for i in range(A.nrows):
            for j in range(A.ncols):
                s = s + A.matrix[i][j] * A.matrix[i][j]
        return root(s)
    #

    #Перегрузка операций для объектов класса матриц:

    #Сложение
    def __add__(self4, other4):
        if (self4.nrows == other4.nrows) and (self4.ncols == other4.ncols):
            
            A = [[0] * self4.ncols for i in range(self4.nrows)]

            for i in range (self4.nrows):
                for j in range(self4.ncols):
                    A[i][j] = self4.matrix[i][j] + other4.matrix[i][j]
        else:
            print("Wrong dimensions! (__add__)")
            A = 0
        return Matrix(A)
    #

    #Вычитание
    def __sub__(self5, other5):
        if (self5.nrows == other5.nrows) and (self5.ncols == other5.ncols):
            
            A = [[0] * self5.ncols for i in range(self5.nrows)]

            for i in range (self5.nrows):
                for j in range(self5.ncols):
                    A[i][j] = self5.matrix[i][j] - other5.matrix[i][j]
        else:
            print("Wrong dimensions! (__sub__)")
            A = 0
        return Matrix(A)
    #

    #Умножение
    def __mul__(self6, other6):
        if (self6.ncols == other6.nrows):
            
            A = [[0] * other6.ncols for i in range(self6.nrows)]

            for i in range (self6.nrows):
                for j in range(other6.ncols):
                    A[i][j] = Matrix.scal_mul(self6.matrix[i], other6.trans()[j])
        else:
            print("Wrong dimensions! (__mul__)")
            A = 0
        return Matrix(A)
    #

    #Вывод:
    def __str__(self7):
        return str(self7.matrix)
    #

    ##


    ###МЕТОД ГАУССА. LUP-РАЗЛОЖЕНИЕ:
    ##
    def Carolus_Fridericus_Gauss_LUP(self, b):

        #
        n = int(self.nrows)
        P = Matrix.eye(self.nrows);
        M = Matrix.eye(n) * Matrix(self.matrix.copy())

        for i in range(self.nrows-1):

            #Поиск опорного элемента в i-ом столбце:
            t =  abs(float(M.matrix[i][i]))
            I = i
            for j in range(i+1, self.nrows):
                if (abs(M.matrix[j][i]) > t):
                    t = abs(float(M.matrix[j][i]))
                    I = j
            M = M.swap(i, I)
            P = P.swap(i, I)

            #Преобразование матрицы M:

            #1)Делим все элементы i-го столбца на опорный элемент:
            for j in range(i+1, n):
                M.matrix[j][i] = float(M.matrix[j][i])/float(M.matrix[i][i])

            #2)Преобразуем все элементы ниже i-ой строки и правее i-го столбца
            for j in range(i+1, self.nrows):
                for k in range(i+1, self.ncols):
                    M.matrix[j][k] = float(M.matrix[j][k]) - float(M.matrix[j][i])*float(M.matrix[i][k])
            #


        #LU-разложение. 
        #M = L - E + U  =>  L + U = M + E
        
        U = Matrix([[0] * self.ncols for i in range(self.nrows)])
        L = Matrix([[0] * self.ncols for i in range(self.nrows)]) + Matrix.eye(self.nrows)
        for i in range(self.nrows):
            for j in range(self.ncols):
                if (j >= i):
                    U.matrix[i][j] = float(M.matrix[i][j])
                else:
                    L.matrix[i][j]= float(M.matrix[i][j])
        #

        #ОБратный ход:
        #1)Ly = Pb
        g = P * b
        y = Matrix(Matrix([[0]*int(self.nrows)]).trans())
        y.matrix[0][0] = g.matrix[0][0]/L.matrix[0][0]
        for i in range(1, self.nrows):
            s = 0
            for j in range(i):
                s = s + L.matrix[i][j]*y.matrix[j][0]
            y.matrix[i][0] = (g.matrix[i][0] - s)/(L.matrix[i][i])
        #

        #2)Ux=y
        x = Matrix(Matrix([[0]*self.nrows]).trans())
        x.matrix[self.nrows-1][0] = y.matrix[self.nrows-1][0] / U.matrix[self.nrows-1][self.nrows-1]

        for i in range(self.nrows-2, -1, -1):
            s = 0
            for j in range(self.nrows - 1, i, -1):
                s = s + U.matrix[i][j]*x.matrix[j][0]
            x.matrix[i][0] = (y.matrix[i][0] - s)/U.matrix[i][i]
        #

        return x
        #

    ##


    ##МЕТОД ХАУСХОЛДЕРА:
    def Alston_Scott_Householder(self, b):

        #Начальные данные:
        n = int(self.nrows)
        m = int(self.ncols)
        Q = Matrix.eye(n)
        R = Matrix(self.matrix.copy())
        #

        #Реализация:
        
        for i in range(n-1):
            
            #Векторы:
            z = Matrix([Matrix.eye(n-i).matrix[0]])
            z = Matrix(Matrix.trans(z))
            y = Matrix([[a[i]] for a in R.matrix[i:n]])
            alph = Matrix.norm_vect(Matrix.trans(y)[0])
            p = Matrix.norm_vect(Matrix.trans(y - Matrix.eye(n-i, alph) * z)[0])
            w = Matrix.eye(n-i, 1/p) * (y - Matrix.eye(n-i, alph) * z)


            #Матрицы:
            Q_0 = Matrix.eye(n-i) - Matrix.eye(n-i, 2) * w * Matrix(Matrix.trans(w))
            A = Matrix.eye(n)
            for k in range (n):
                for j in range(n):
                    if (k >= n - Q_0.nrows) and (j >= n - Q_0.nrows):
                        A.matrix[k][j] = float(Q_0.matrix[k - (n - Q_0.nrows)][j - (n - Q_0.nrows)])
            Q_0 = A
            R = Q_0 * R
            Q = Q * Q_0

        
        #y = Q'*b
        Y = Matrix(Matrix.trans(Q)) * b

        #R*x = y 
        x = Matrix(Matrix.trans(Matrix([[0]*n])))
        x.matrix[n-1][0] = Y.matrix[n-1][0] / R.matrix[n-1][n-1]

        for i in range(n-2, -1, -1):
            s = 0
            for j in range(n - 1, i, -1):
                s = s + R.matrix[i][j]*x.matrix[j][0]
            x.matrix[i][0] = (Y.matrix[i][0] - s)/R.matrix[i][i]
        #
        #
        print("Результат в методе Хаусхолдера: ")
        return x
        #
    ##




        ###МЕТОД ПРОСТЫХ ИТЕРАЦИЙ:
    def MPI(self, d, eps = 10**(-3)):

        #Начальные данные
        n = int(self.nrows)
        A = Matrix.eye(n) * Matrix(self.matrix.copy())
        mu = 1 / Matrix.minim([Matrix.Norma_1(A), Matrix.Norma_2(A), Matrix.Norma_3(A)])
        B = Matrix.eye(n) - Matrix.eye(n, mu) * A
        #

        #Проверка условий сходимости и её обеспечение путём приведения системы домножением слева на A':
        norm = Matrix.minim([Matrix.Norma_1(B), Matrix.Norma_2(B), Matrix.Norma_3(B)])
        if (norm >= 1):
            A_ = Matrix(Matrix.trans(A))
            A = A_ * A
            d = A_ * d
            mu = 1 / Matrix.minim([Matrix.Norma_1(A), Matrix.Norma_2(A), Matrix.Norma_3(A)])
            B = Matrix.eye(n) - Matrix.eye(n, mu) * A
        #

        #Считаем вектор c:
        c = Matrix.eye(n, mu) * d
        #
        k = 0
        #Ещё одна проверка сходимости нормы
        norm = Matrix.minim([Matrix.Norma_1(B), Matrix.Norma_2(B), Matrix.Norma_3(B)])
        if (norm < 1):
            x0 = Matrix(c.matrix.copy())
            x = x0 + Matrix([[eps * (norm / (1 - norm))]] * n)
            while (Matrix.norm_vect(Matrix.trans(x - x0)[0]) * (norm / (1 - norm)) >= eps) and (k < 1000):
                x0 = Matrix(x.matrix.copy())
                x = B * x + c
                k +=1
                if (k == 1000):
                    print("The number of iterations in MPI method has reached ten thousand")
        else:
            x = Matrix(c.matrix.copy())
            while (Matrix.norm_vect(Matrix.trans(A * x - d)[0]) >= eps) and (k <1000):
                x = B * x + c
                k +=1
                if (k == 1000):
                    print("The number of iterations in MPI method has reached one thousand")
        print("Количество итераций: ", k)
        print("Результат в Методе Простых Итераций:")
        return x
        #
    ###

###МЕТОД ЗЕЙДЕЛЯ
    def Philipp_Ludwig_von_Seidel(self, b, eps = 10**(-3)):

        #Начальные данные
        n = int(self.nrows)
        A = Matrix.eye(n) * Matrix(self.matrix.copy())
        for i in range(n):
            if (A.matrix[i][i] == 0):
                for j in range(n):
                    if ((int(A.matrix[j][i]) != 0) and (int(A.matrix[i][j]) != 0) and (int(A.matrix[i][i]) == 0)):
                        A = A.swap(i, j)
                        b = b.swap(i, j)
        C = Matrix([[0] * n for i in range(n)])
        d = Matrix([[0] for i in range(n)])
        #

        #Построение C и d:
        for i in range(n):
            for j in range(n):
                if (i != j):
                    C.matrix[i][j] = float(-A.matrix[i][j]/A.matrix[i][i])
            d.matrix[i] = [b.matrix[i][0]/A.matrix[i][i]].copy()
        #

        #Проверка условий сходимости и её обеспечение путём приведения системы домножением слева на A':
        norm = Matrix.minim([Matrix.Norma_1(C), Matrix.Norma_2(C), Matrix.Norma_3(C)])
        if (norm >= 1):
            A_ = Matrix(Matrix.trans(A))
            A = A_ * A
            b = A_ * b
            #Построение C и d снова:
            for i in range(n):
                for j in range(n):
                    if (i != j):
                        C.matrix[i][j] = float(-A.matrix[i][j]/A.matrix[i][i])
                d.matrix[i] = [b.matrix[i][0]/A.matrix[i][i]].copy()
        #
        #


        #Ещё одна проверка сходимости нормы
        x =  Matrix.eye(n) * Matrix(d.matrix.copy())
        k = 0
        while (Matrix.norm_vect(Matrix.trans(A * x - b)[0]) >= eps) and (k < 1000):
            for i in range(n):
                x.matrix[i][0] = Matrix.scal_mul(C.matrix[i], Matrix.trans(x)[0]) + float(d.matrix[i][0])
            k += 1
            if (k == 1000):
                print("The number of iterations in MPI method has reached ten thousand")

        print("Количество итераций: ", k)
        print("Результат в методе Зейделя:")
        return x
        #
    ###

#######Нормальные полиномы:
def MNK_norm(x, X, F, n):
    E = [] 
    for i in range(len(X)):
        E.append([X[i]**j for j in range(n+1)])
    F = [[F[j]] for j in range(len(F))]
    E = Matrix(E)
    F = Matrix(F)

    a = Matrix.Carolus_Fridericus_Gauss_LUP(Matrix(E.trans())*E, Matrix(E.trans())*F)

    return sum([a.matrix[i][0]*(x**i) for i in range(len(a.matrix))])


#Реккурентная функция ортогональных полиномов
def q_j(x, X, j):
    if (j > 1):
        a = ( sum( [ X[i] * ((q_j(X[i], X, j-1))**2) for i in range(len(X))] ) / sum( [ ((q_j(X[i], X, j-1))**2) for i in range(len(X))] ) )
        b = ( sum( [ X[i] * (q_j(X[i], X, j-1)) * (q_j(X[i], X, j-2)) for i in range(len(X))] ) / sum( [ ((q_j(X[i], X, j-2))**2) for i in range(len(X))] ) )
        return x*q_j(x, X, j-1) - a * q_j(x, X, j-1) - b * q_j(x, X, j - 2)
    elif (j == 0):
        return 1
    elif (j == 1):
        return x - sum(X)/(len(X))

#Ортогональные полиномы
def MNK_ort(x, X, F, n):
    a = []
    for i in range(n+1):
        a.append(sum( [ q_j(X[j], X, i)*F[j] for j in range(len(X))] ) / sum( [ q_j(X[t], X, i)**2 for t in range(len(X))] ))
        print('Iteration', i, 'from', n)
    return sum( [ a[i]*q_j(x, X, i) for i in range(len(a)) ] )


#Ортогональные полиномы:
def MNK_ort_(x, X, F, n):
    Q = []
    Q.append([1]*len(X))
    Q.append([t - sum(X)/(len(X)) for t in X])
    for i in range(2, n+1):
        a = ( sum( [ X[j] * (Q[i-1][j])**2 for j in range(len(X))] ) / sum( [ ((Q[i-1][j])**2) for j in range(len(X))] ) )
        b = ( sum( [ X[j] * (Q[i-1][j]) * (Q[i-2][j]) for j in range(len(X))] ) / sum( [ ((Q[i-2][j])**2) for j in range(len(X))] ) )
        Q.append([X[j]*Q[i-1][j] - a * Q[i-1][j] - b * Q[i-2][j] for j in range(len(X))])
        #
    A = []
    for i in range(n+1):
        A.append(sum( [ Q[i][j]*F[j] for j in range(len(X))] ) / sum( [ Q[i][t]**2 for t in range(len(X))] ))
    q = [[1]*len(x)]
    q.append([t - sum(X)/(len(X)) for t in x])
    for i in range(2, n+1):
        a = ( sum( [ X[j] * (Q[i-1][j])**2 for j in range(len(X))] ) / sum( [ ((Q[i-1][j])**2) for j in range(len(X))] ) )
        b = ( sum( [ X[j] * (Q[i-1][j]) * (Q[i-2][j]) for j in range(len(X))] ) / sum( [ ((Q[i-2][j])**2) for j in range(len(X))] ) )
        q.append([x[j]*q[i-1][j] - a * q[i-1][j] - b * q[i-2][j] for j in range(len(x))])
    return [sum([ A[j]*q[j][i] for j in range(len(A)) ]) for i in range(len(x))], sum([ (sum( [ A[j]*Q[j][i] for j in range(len(A)) ]) - F[i])**2 for i in range(len(X)) ])



                        ######################################################################
N = 2                                   #Порядок полинома - 1
N_1 = 10                                #Количество различных точек
A = -2                                  #Границы
B = 3                                   #
i = 5                                   #Количество повторяющихся 
eps = 0.7                               #Шумовые отклонения
f = lambda x: np.sin(x)

x = np.linspace(A, B, N_1)
X = []
for j in range(len(x)):
    X += [x[j]]*i

F = []
for j in range(len(x)):
    F.append(f(x[j]))
    for k in range(i-1):
        F.append(2*eps*(np.random.rand()-0.5) + f(x[j]))


x = np.linspace(A, B, 100)

fig, axs = plt.subplots(5)
for i in range(5):
    axs[i].plot(x, f(x), label = 'f(x)')
    axs[i].plot(x, MNK_ort_(x, X, F, i+1)[0], label = 'n = '+str(i+1))
    for j in range(len(X)):
        axs[i].scatter(X[j], F[j], marker = 'o', s = 20)
    axs[i].legend()
    axs[i].grid()
fig.suptitle('Ортогональные полиномы')

#df = []
#for j in range(5):
#    df.append([j+1, sum( [ (MNK_norm(X[i], X, F, j+1) - F[i])**2 for i in range(len(X))] ), MNK_ort_(x, X, F, j+1)[1]])
#data = pd.DataFrame(df, index = [j+1 for j in range(5)], columns = ['n', 'Сумма квадратов ошибок для нормальных уравнений', 'Сумма квадратов ошибок для ортогональных полиномов'])
#data.to_excel("output.xlsx")  
plt.show()