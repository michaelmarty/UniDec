from TailFragments2 import *
import textwrap
# print(predict_tail_fragments("14:1;2O"))
print(textwrap.fill(str(predict_tail_fragments("18:0;2O")), width=100))
# print(predict_tail_fragments("O-14:0", hgfrags=["224.07"]))
# print(predict_tail_fragments("8:0;3O", hgfrags=["179.06"], hgnls=[46.01,208.06]))
# print(predict_tail_fragments("12:1;(2OH)"))