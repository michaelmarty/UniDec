# noinspection PyUnresolvedReferences
import sys
from unidec.engine import *


def main(*args, **kwargs):
    print(sys.argv)
    print("Running Command Line UniDec")
    UniDec(*sys.argv)


if __name__ == '__main__':
    main()

