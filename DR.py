# -*- coding: utf8 -*-
from numpy import *
import numpy as np
from termcolor import colored


EPSILON_P = 0.000000000001  # точность округления
EPSILON_M = -0.000000000001  # точность округления

def color_pick(_i):
    """Выбираем цвет"""
    color = ['red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'grey']
    index = _i - int(_i / 8) * 8
    return color[index]


def read_data():
    """Чтение данных их файла"""

    file_name = input("Введите номер примера: ")

    # открытие файла
    with open(file_name + ".dat") as ifs:
        lines = ifs.readlines()

    # кол. уравнений  и кол. переменных
    lines[0].split()
    m = int(lines[0][0])
    n = int(lines[0][2])

    # создание пустых списков
    a = np.zeros((m, n + 2), dtype=float)
    s = np.zeros((m, 1), dtype=str)

    # считывание и формирование СЛАУ
    for i in range(m):
        lines[i] = lines[i + 1].split()

        for j in range(n + 2):
            if j < n:
                a[i, j] = float(lines[i][j])
            elif j == n:
                s[i, 0] = lines[i][j]
            else:
                a[i, n + 1] = float(lines[i][j])

    return a, s


def print_list(_list, _flag):
    """Вывод списка"""

    if _flag == 0:
        print(colored("\nРасширеная матрица: ", color_pick(7)))
        print(colored(_list, color_pick(7)))
        print("")
    elif _flag == 1:
        print(colored("\nРасширенная матрица с ИБ: ", color_pick(7)))
        print(colored(_list, color_pick(7)))
        print("")
    elif _flag == 2:
        print(colored("Полученная симплекс таблица: ", color_pick(7)))
        print(colored(_list, color_pick(7)))
        print("")


def print_vector(_vector, _flag):
    """Вывод вектора"""

    if _flag == 0:
        print(colored("Базисные переменные: ", color_pick(7)))
        print("\t".join([colored(str(k), color_pick(0)) for k in _vector]))
        print("")
    elif _flag == 1:
        print(colored("Свободные переменные: ", color_pick(7)))
        print("\t".join([colored(str(k), color_pick(1)) for k in _vector]))
        print("")
    elif _flag == 2:
        print(colored("Первое базисное решение: ", color_pick(7)))
        print("[", end='')
        print("\t".join([colored(str(k), color_pick(3)) for k in _vector]), end='')
        print("]")
    elif _flag == 3:
        print(colored("Полученное базисное решение: ", color_pick(7)))
        print("[", end='')
        print("\t".join([colored(str(k), color_pick(3)) for k in _vector]), end='')
        print("]")
    elif _flag == 4:
        print(colored("Функция цели: ", color_pick(7)))
        print("[", end='')
        print("\t".join([colored(str(k), color_pick(3)) for k in _vector]), end='')
        print("]")


def print_solve(_i, _j):
    """Вывод решения по индексам"""
    if _i == -1 and _j == -1:
        print(colored("Нет допустимого базисного решения системы!", color_pick(0)))
        return False
    elif _i == -2 and _j == -2:
        print(colored("Область допустимых решений не ограничена!", color_pick(0)))
        return  False
    elif _i == -3 and _j == -3:
        print(colored("Область допустимых решений не ограничена!", color_pick(0)))
        return False
    else:
        print(colored("Индексы разрешающего элемента: [i=%d, j=%d]\n" % (_i + 1, _j + 1) , color_pick(7)))
        return True


def to_equations(_a, _s):
    """Вводим новые переменные x[n + 1...m]"""
    loc = len(_a[0]) - 2

    for i in range(len(_a)):
        # вводим +x
        if _s[i, 0] == "<":
            _a[i, loc] = 1
        # вводим -x
        elif _s[i, 0] == ">":
            _a[i, loc] = -1

    for i in range(len(_a)):
        if _a[i, len(_a[0]) - 1] < 0:
            _a[i] *= -1

    return _a


def vector_base(_m, _n):
    """Создание вектора базисных"""

    vector = []
    vector_ind = []
    for i in range(_n, _m + _n):
        vector.append("x" + str(i + 1))
        vector_ind.append(i)

    return vector, vector_ind


def vector_empty(_n):
    """Создание вектора свободных"""

    vector = []
    vector_ind = []

    for i in range(_n):
        vector.append("x" + str(i + 1))
        vector_ind.append(i)

    return vector, vector_ind


def vector_base_a(_m, _n ,_an, _i_ab):
    """Создание вектора базисных с ИБ"""
    vector = []
    vector_ind = []

    for i in range(_m):
        if _i_ab[i] == 1:
            vector.append("x" + str(_an - 1 + i))
            vector_ind.append(_an - 2 + i)
        else:
             vector.append("x" + str(_n + i + 1))
             vector_ind.append(_n + i)

    return vector, vector_ind


def vector_empty_a(_n, _m, _i_ab):
    """Создание вектора свободных c ИБ"""
    vector = []
    vector_ind = []

    for i in range(_n):
         vector.append("x" + str(i + 1))
         vector_ind.append(i)
    for i in range(_n, _n + _m):
        if _i_ab[i - _n] != 0:
            vector.append("x" + str(i + 1))
            vector_ind.append(i)

    return vector, vector_ind


def first_plan(_a):
    """Создание опорного плана"""

    plan = []
    n = len(_a[0]) - 2
    k = n

    for i in range(len(_a) + n):
        if i < k:
            plan.append(0)
        else:
            plan.append(_a[i - k, n + 1] * _a[i - k, n])

    return plan


def new_base_plan(_a, _bxi, _exi):
    """Создание опорного плана"""
    plan = []

    for i in range(len(_exi) + len(_bxi)):
        plan.append(0)

    for i in range(len(_bxi)):
        plan[_bxi[i]] = _a[i][0]

    return plan


def test_plan(_p):
    """Проверка базисного плана на допустимость"""

    for i in range(len(_p)):
        if _p[i] < 0:
            print(colored("Базисные переменные отрицательный!", color_pick(4)))
            return False

    print(colored("Получено допустимое базисное решение.", color_pick(4)))
    return True


def add_artificial_base(_a):
    """Дополнение искусственым базисом"""

    m = len(_a)
    n = len(_a[0]) - 2
    a = np.zeros((m, n + 3), dtype=float)

    for i in range(m):
        for j in range(n + 1):
            a[i, j] = _a[i, j]
            a[i, len(a[0]) - 1] = _a[i, len(_a[0]) - 1]

            if _a[i, len(_a) - 1] < 0 or _a[i, len(_a)] < 0:
                a[i, n + 1] = 1
            else:
                a[i, n + 1] = 0

    return a


def calc_n(_a, _m, _n):
    """Подсчет количества переменных с базисными"""
    n = _m + _n
    ab = 0
    ind_ab = []

    for i in range(len(_a)):
        if _a[i, _n + 1] == 1:
            ab += 1
            ind_ab.append(1)
        else:
            ind_ab.append(0)

    return n + ab, ab, ind_ab


def new_f(_a):
    """Подсчет новой функции цели"""
    return _a[len(_a) - 1]


def swap_x(_bx, _ex, _bxi, _exi, _ri, _rj, _t_iab):
    """смена x"""

    tmp, tmp_ind = _bx[_ri], _bxi[_ri]
    _bx[_ri], _bxi[_ri] = _ex[_rj - 1], _exi[_rj - 1]
    _ex[_rj - 1], _exi[_rj - 1] = tmp, tmp_ind

    ex = []
    exi = []

    if _t_iab[_ri] == 0:
        ex = _ex
        exi = _exi
    else:
        for i in range(len(_ex)):
            if i != _ri:
                ex.append(_ex[i])
                exi.append(_exi[i])

    _t_iab[_ri] = 0

    return _bx, ex, _bxi, exi, _t_iab


def delete_pill(_a, _i):
    """Удаление столбца"""
    return np.delete(_a, _i, axis=1)


def simplex_table(_a, _n):
    """Создание симплекс таблицы для смены базиса"""

    m = len(_a)
    n = len(_a[0])
    emp = 0

    # подсчет свободных доп. x
    for i in range(m):
        if _a[i, n - 2] == 1:
            emp += 1

    new_n = _n + emp + 1

    table = np.zeros((m + 1, new_n ), dtype=float)
    fg = np.zeros((emp, new_n), dtype=float)

    # создание начальной таблицы
    for i in range(m):
        table[i, 0] = _a[i, n - 1]
        for j in range(1, _n + 1):
            table[i, j] = _a[i, j - 1]

    # создание масива ei
    num = 0
    ind_str = []
    for i in range(m):
        if _a[i, n - 2] == 1:
            ind_str.append(i)
            for j in range(new_n):
                if j <= _n:
                    fg[num, j] = table[i, j]
                else:
                    fg[num, j + num] = -1
                    break
            num += 1

    # ввод базиса в таблицу
    for i in range(emp):
        table[ind_str[i]] = fg[i]

    # подсчет функции цели
    f = np.zeros((1, new_n), dtype=float)

    for i in range(emp):
        f += fg[i]* -1
    table[len(table) - 1] = f

    return table, f


def new_simplex_table(_a, _ri, _rj, _i_ab, _t_iab, _au_iab):
    """Создание новой симлекс таблицы"""
    m = len(_a)
    n = len(_a[0])

    new_table = np.zeros((m, n), dtype=float)

    re = 1.0 / _a[_ri, _rj]
    new_table[_ri, _rj] = re

    # подсчет строки
    for j in range(n):
        if j != _rj:
            new_table[_ri, j] = _a[_ri, j] * re

    # подсчет столбца
    for i in range(m):
        if i != _ri:
            new_table[i, _rj] = (_a[i, _rj] * re) * -1.0

    # подсчет элементов
    for i in range(m):
        for j in range(n):
            if j != _rj and i != _ri:
                new_table[i, j] = _a[i, j] - ((_a[_ri, j] * _a[i, _rj]) * re)

    # удаление ИБ свободной
    if _t_iab[_ri] == 0 and _i_ab[_ri] == 1 and _au_iab[_ri] == 0:
        new_table = delete_pill(new_table, _rj)
        _au_iab[_ri] = 1

    # зануление малых величин
    for i in range(len(new_table)):
        for j in range(len(new_table[0])):
            if 0 < abs(new_table[i, j]) < EPSILON_P:
                new_table[i, j] = 0

    f = np.zeros((1, len(new_table[0])), dtype=float)
    for j in range(len(new_table[0])):
        f[0, j] = new_table[len(new_table) - 1, j]

    return new_table, f, _au_iab


def find_re(_a):
    """Поиск индексов разрешающего элемента"""
    ri, rj = -1, -1

    m = len(_a)
    n = len(_a[0])
    min_m = []
    f_str = []

    for i in range(1, n):
        f_str.append(_a[m - 1, i])

    while True:
        mini = max(f_str)
        for i in range(len(f_str)):
            if f_str[i] < 0:
                if f_str[i] < mini:
                    mini = f_str[i]
                    rj = i + 1

        if rj == -1:
            return -1, -1

        ind_min = []
        for i in range(m):
            if _a[i, rj] > 0:
                min_m.append(_a[i, 0] / _a[i, rj])
                ind_min.append(i)

        if len(min_m) == 0 and len(f_str) != 0:
            f_str.remove(f_str[rj])
            return -3, -3
            #continue
        elif len(min_m) == 0 and len(f_str) == 0:
            return -2, -2
        else:
            ri = ind_min[min_m.index(min(min_m))]
            return ri, rj

def main():
    """1 часть, если есть первый опорный план"""
    a, signs = read_data() # чтение неравенств
    a = to_equations(a, signs)  # дополнение до уравнений

    m = len(a)
    n = len(a[0]) - 2

    # формирование строковых списков
    empty_x, empty_x_ind = vector_empty(n)
    base_x, base_x_ind = vector_base(m, n)

    # вывод матрицы и векторов
    print_list(a, 0)
    print_vector(base_x, 0)
    print_vector(empty_x, 1)

    # получение 1ого о.п
    base_plan = first_plan(a)
    print_vector(base_plan, 2)

    if test_plan(base_plan):
        return 0

    print("".join([colored("-", color_pick(5)) for i in range(100)]), end='')

    """2 часть, метод ИБ"""
    a = add_artificial_base(a)

    an, ab, i_ab = calc_n(a, m, n)

    tempura = []
    au_iab = []
    for i in range(len(i_ab)):
        tempura.append(i_ab[i])
        au_iab.append(0)

    # формирование строковых списков
    empty_x, empty_x_ind = vector_empty_a(n, m, tempura)
    base_x, base_x_ind = vector_base_a(m, n, an, tempura)

    # вывод матрицы и векторов
    print_list(a, 1)
    print_vector(base_x, 0)
    print_vector(empty_x, 1)

    # создание симлекс таблицы
    a, f = simplex_table(a, n)
    print_list(a, 2)
    print_vector(f, 4)

    index_i, index_j = find_re(a)

    if not print_solve(index_i, index_j):
        return 0

    base_x, empty_x, base_x_ind, empty_x_ind, tempura = swap_x(base_x, empty_x, base_x_ind, empty_x_ind, index_i, index_j, tempura)


    count = 1
    while True:
        a, f, au_iab = new_simplex_table(a, index_i, index_j, i_ab, tempura, au_iab)

        print("Итерация: " + str(count))
        print_list(a, 2)
        print_vector(f, 4)

        index_i, index_j = find_re(a)

        if index_i == -1 and index_j == -1:

            mm = f[0, 0]
            for k in range(len(f[0])):
                if f[0, k] <= mm:
                    mm = f[0, k]

            # Если нашли решение
            if mm >= 0 and f[0, 0] == 0:
                print(colored("\nНайдено допустимое базисное решение!", color_pick(4)))
                print_vector(base_x, 0)
                print_vector(empty_x, 1)
                # Построение базисного плана
                base_plan = new_base_plan(a, base_x_ind, empty_x_ind)
                print_vector(base_plan, 3)
                return 0
            elif f[0, 0] < 0:
                print(colored("\nИсходная задача не имеет решения(система условий этой задачи несовместна)!", color_pick(4)))
                return 0
        else:
            if print_solve(index_i, index_j):
                base_x, empty_x, base_x_ind, empty_x_ind, tempura = swap_x(base_x, empty_x, base_x_ind, empty_x_ind, index_i, index_j, tempura)
                count += 1
            else:
                return 0



if __name__ == '__main__':
    main()
