# noinspection PyUnresolvedReferences
import sys
from unidec.engine import *


def main(*args, **kwargs):
    print("\nRunning Command Line UniDec")
    print("Args:", sys.argv)
    UniDec(*sys.argv)


if __name__ == '__main__':
    main()

