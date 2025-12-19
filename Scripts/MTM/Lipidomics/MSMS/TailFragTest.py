from TailFragmentsAssigner import *
import textwrap
# print(predict_tail_fragments("14:1;2O"))
frags, nls = predict_tail_fragments("18:1;(2OH)", classname="", mode="positive")
print(textwrap.fill(str(frags), width=100))
print("\n")
print(textwrap.fill(str(nls), width=100))
# print(predict_tail_fragments("O-14:0", hgfrags=["224.07"]))
# print(predict_tail_fragments("8:0;3O", hgfrags=["179.06"], hgnls=[46.01,208.06]))
# print(predict_tail_fragments("12:1;(2OH)"))
input()