import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style(style='white')

def main():

    input_dict = {'fit_range': [2, 8]}

    min_fit, max_fit = input_dict['fit_range']

    X = my_array([1,3,4,5,6])
    X.randomize()
    print(X)

    X = my_array(list(np.linspace(0,10, 100)))
    Y = my_array(list(np.linspace(1,11,100)))

    X.randomize()
    X.sort()
    Y.randomize()

    # plot_x_vs_y(X,Y)

    X_restricted, Y_restricted = X.restrict(Y, min_fit, max_fit)

    print('classes X: {} Y: {}'.format(type(X), type(Y)))

    plot_x_vs_y(X_restricted, Y_restricted)

    a, b = compute_a_b(X_restricted, Y_restricted)  # coef of line

    print(a,b)

    x0 = 0
    y0 = b

    plot_extrapolation(X,Y, X_restricted, Y_restricted)


class my_array(np.ndarray):

    # def __new__(cls):
    #     cls = super.__new__(cls=np.ndarray)

    def __new__(cls, content):
        shape = np.shape(content)
        return super().__new__(cls, shape)

    def __init__(self, content):
        self=np.array(content)

    def randomize(self):
        len_s  = len(self)
        for i in range(len_s):
            self[i] = self[i] + np.random.random()
#        print('this is self at the exit of the function', self)

    def restrict(self, Y, a, b):
        target_i = np.where(np.logical_and(self> a, self < b))[0]
        selected_r = self[target_i]  # selected interval to compute epsilon
        selected_energies = Y[target_i]  # TODO: all energies to IP/EA
        return selected_r, selected_energies

def plot_x_vs_y(x, y, label_x='x', label_y='y', LineStyle='solid', marker='.'):
    fig, ax = plt.subplots()
    ax.plot(x,y,LineStyle=LineStyle, marker=marker)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    # plt.show()

def plot_extrapolation(x, y, x_, y_, label_x='x', label_y='y', LineStyle='None', marker='.', xmin = 0, ymin = 0):
    a,b = compute_a_b(x_, y_)
    x0, y0 = 0, b
    x1, y1 = max(x), a*max(x) + b
    fig, ax = plt.subplots()
    ax.plot(x,y, LineStyle=LineStyle ,marker=marker)
    ax.plot(x_, y_, LineStyle=LineStyle, marker='o')
    ax.plot([x0, x1] , [y0, y1])
    ax.plot([x0] , [y0], marker='x')
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    ax.grid()
    #
    # sns.lineplot(x=x_, y=y_, err_style="bars", ci=68)
    #
    ax.set_xlim(left=xmin)
    ax.set_ylim(bottom=ymin)
    plt.show()

def compute_a_b(X,Y):
    a,b = np.polyfit(X,Y,1)
    return a,b

def compute_a_b(X,Y):
    a,b = np.polyfit(X,Y,1)
    return a,b

if __name__ == '__main__':
    main()