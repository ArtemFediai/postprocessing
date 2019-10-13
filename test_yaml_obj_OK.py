"""
test ruamel yaml
"""

from ruamel.yaml import YAML
import numpy as np
from ruamel.yaml import load


def main():
    # make a library

    a = np.array([1.0, 1.244840490])
    b = np.array([2.0, 4.56654])
    c = {0: a.tolist(), 1: b.tolist()}

    # make it an object instance
    my_instance = MyClass(c)

    # register the class in yaml
    yaml = YAML()
    yaml.register_class(MyClass)

    # get fid, dump the class
    fid = open("test.yml", "w")
    yaml.dump(my_instance, fid)
    fid.close()

    # load just created class
    fid = open("test.yml", "r")
    yaml1 = YAML(typ="safe", pure=True)  # ini
    yaml1.register_class(MyClass)
    read_out_obj = yaml1.load(fid)
    fid.close()

    # check if classes attributes are the same
    AreTheyTheSame = (my_instance.__dict__ == read_out_obj.__dict__)

    pass


class MyClass:
    def __init__(self, c):
        self.c = c


if __name__ == '__main__':
    main()
