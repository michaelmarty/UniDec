seq = "000100110101111"
s = list(seq)
ints = [eval(i) for i in s]
inv_seq = [abs(i-1) for i in ints]
inv_seq = ''.join(map(str, inv_seq))
print(inv_seq)



#seq = seq.split()
#inv_seq = abs(seq - 1)
#print(inv_seq)