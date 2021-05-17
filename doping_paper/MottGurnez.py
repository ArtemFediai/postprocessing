import os.path
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.linear_model


def main():
    print("this is Mott Gurnetz low")
    my_path = '/home/artem/Desktop/295.csv'
    # <--  change it to whatever you want. It must have J, V columns in CI (A/m^2 and V)

    path_to_csv = os.path.abspath(my_path)
    print(path_to_csv)

    J, V = [], []

    # --> just another method to import values
    # with open(file=path_to_csv, mode='r') as stream:
    #     content = csv.reader(stream, delimiter=',')
    #     for i, line in enumerate(content):
    #         V.append(float(line[0]))
    #         J.append(float(line[1]))
    #         print(line)
    # print('current', J)
    # print('voltage', V)
    # <-- just another ...

    pd_data = pd.read_csv(filepath_or_buffer=path_to_csv, header=None)
    new_pd = pd_data.sort_values(by=1)
    print(pd_data)
    print(new_pd)

    plt.Figure()
    # fig, ax = plt.subplots()
    # ax.plot(V, J, label='normal import')
    new_pd.plot(x=0, y=1, logx=True, logy=True, title='aNPD', grid=True, legend=True, label='aNPD')
    new_pd.columns = ['V', 'J']  # give names to J and V columns
    new_pd['V^2'] = [x*x for x in new_pd['V']]  # create and append V^2
    print(new_pd)

    new_pd.plot(x='V^2', y='J', linestyle=None, marker='.')
    X, Y = np.array(new_pd['V^2']), np.array(new_pd['J'])
    reg = sklearn.linear_model.LinearRegression().fit(X.reshape(-1, 1), Y)
    y_predict = reg.predict(X.reshape(-1,1))
    plt.plot(X, y_predict, label='pred')
    plt.legend()

    print(y_predict)
    alpha = reg.coef_  # slope of the J vs V^2
    print(alpha)
    coef = mg_low_coef()  # coef in the formula Mott Gurnetz
    mu = alpha/coef  #  this follows from the Mott Gurnetz formula
    R2 = reg.score(X.reshape(-1, 1), Y)
    print(f'determination factor R^2: {R2}')
    print(f'mobility: {mu}')
    plt.show()


def mg_low_coef(eps_r=3, l=200*1E-9):
    """
    J = 9/8 * eps_r * eps_0 * mu * V**2 * L**-3
    :param eps_r: relative dielectric permittivity
    :param d: thickness of the layer
    :param EPS0: known universal constant
    :return: coefficient in the formula
    """
    EPS0 = 8.854E-12
    coef = 9/8 * eps_r * EPS0 * l**(-3)
    return coef


if __name__ == '__main__':
    main()
