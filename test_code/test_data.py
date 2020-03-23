import numpy as np
import pandas
df = pandas.read_csv('PhiNN.txt',header=None)
a = df.to_numpy()
# a =  np.char.replace(a,'i','j').astype(np.complex) 

b = []
for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        b.append(complex(a[i,j])) 
c = np.asarray(b)
c = np.reshape(c, (2,2))

w, v = np.linalg.eig(c)

for i in range(2):
    print("\n")
    print(v[:,i])
    # # print(i+1)
    # # print(w[i])
    # # print("\n")
    # # print(v[i])
    # print(w[i]*v[i]-c.dot(v[i]))
