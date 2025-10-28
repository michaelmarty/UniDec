import sys

sys.path.append("C:\\Python\\UniDec")
import sklearn.cluster as skl

fname = "GroEL_td.txt"

newdir = os.path.join("C:\\Python\\UniDec", "TestSpectra")

data = np.loadtxt(os.path.join(newdir, fname))
data[:, 1] /= np.amax(data[:, 1])
#data=limit_data(data,0,4000)

print(auto_peak_width(data))

cs=np.transpose([data[:,0],np.cumsum(data[:,1])])
print(cs.shape)
cs[:,1]/=np.amax(cs[:,1])
cs[:,0]/=np.amax(cs[:,0])

skl.DBSCAN()

plt.plot(cs[:,0],cs[:,1])
#plt.show()

input = cs[:,0].reshape(-1,1) #np.transpose([data[:, 1], cwt])
print("clustering")
#db = skl.DBSCAN(eps=0.01, min_samples=10).fit(input)
db = skl.KMeans().fit(input)
print("done")
labels = db.labels_
print(np.unique(labels))


core_samples_mask = np.zeros_like(db.labels_, dtype=bool)


unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    xy = input[class_member_mask ]
    x=cs[class_member_mask ,0]
    # plt.plot(xy[:, 0], xy[:, 1], 'o', color=col, markersize=14)
    plt.plot(x, xy, 'o', color=col, markersize=6)


print(np.unique(labels))
plt.show()

exit()


data1=limit_data(data,11160,11200)
data2=limit_data(data,1280,1290)
data3=limit_data(data,11200,11300)
data4=limit_data(data,4000,4200)


plt.figure()
plt.subplot(411)
plt.plot(data1[:,0],data1[:,1])
plt.subplot(412)
plt.plot(data2[:,0],data2[:,1])
plt.subplot(413)
plt.plot(data3[:,0],data3[:,1])
plt.subplot(414)
plt.plot(data4[:,0],data4[:,1])



data1=autocorr(data1)[0]
data2=autocorr(data2)[0]
data3=autocorr(data3)[0]
data4=autocorr(data4)[0]

'''
data1=fft(data1)
data2=fft(data2)
data3=fft(data3)
data4=fft(data4)
'''

plt.figure()
plt.subplot(411)
plt.plot(data1[:,0],data1[:,1])
plt.subplot(412)
plt.plot(data2[:,0],data2[:,1])
plt.subplot(413)
plt.plot(data3[:,0],data3[:,1])
plt.subplot(414)
plt.plot(data4[:,0],data4[:,1])
plt.show()
exit()


dat = win_fft_grid_single(data, .1, 1, 5, diffrange=[0.01,20])

zgrid = dat[:, 2]
xvals = np.unique(dat[:, 0])
yvals = np.unique(dat[:, 1])
xlen = len(xvals)
ylen = len(yvals)
newgrid = np.reshape(zgrid, (xlen, ylen))
print(newgrid.shape, len(xvals), len(yvals))
plt.figure(figsize=(18,12))
plt.subplot(211)
plt.contourf(xvals, yvals, newgrid.transpose())
plt.subplot(212)
plt.plot(data[:,0],data[:,1])
plt.xlim(np.amin(xvals),np.amax(xvals))
plt.show()

exit()

b1 = data[:, 0] > 1000
data = data[b1]

kernel = conv_peak_shape_kernel(data[:, 0], 0, 15)
cdata = cconv(data[:, 1], kernel)
data2 = np.transpose([data[:, 0], cdata])
data2[:, 1] /= np.amax(data2[:, 1])

cwt, wavelet = single_cwt(data[:, 1], 6, wavelet_type="1DG")
plt.plot(wavelet)
plt.show()

out = fft(data2)
print(out.shape, cwt.shape)

r1 = safedivide(data2[:, 1], data[:, 1])
sr1 = np.sort(r1)
sd = np.sort(data[:, 1])
sd2 = np.sort(data2[:, 1])
# plt.plot(1./data2[:,0],data2[:,1])


input = np.transpose([data[:, 1], cwt])
print("clustering")
db = skl.DBSCAN(eps=0.001, min_samples=100).fit(input)
print("done")
labels = db.labels_
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True

unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    xy = input[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', color=col, markersize=14)

    xy = input[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', color=col, markersize=6)

print(labels)
plt.show()
# plt.figure()
'''
plt.plot(sr1)
plt.plot(sd)
plt.plot(sd2)'''
# plt.scatter(data[:,1],cwt)
# plt.xlim(5500,8000)

b2 = labels == 0

plt.figure()
plt.plot(data[b2, 0], data[b2, 1])
plt.plot(data[:, 0], data[:, 1], linestyle="--")
plt.plot(data[:, 0], cwt, ls=":")
plt.show()

'''

if True:
    import numpy as np
    import time
    import theano
    A = np.random.rand(10000,10000).astype(theano.config.floatX)
    B = np.random.rand(10000,10000).astype(theano.config.floatX)
    np_start = time.time()
    AB = A.dot(B)
    np_end = time.time()
    X,Y = theano.tensor.matrices('XY')
    mf = theano.function([X,Y],X.dot(Y))
    t_start = time.time()
    tAB = mf(A,B)
    t_end = time.time()
    print("NP time: %f[s], theano time: %f[s] (times should be close when run on CPU!)" %(
                                               np_end-np_start, t_end-t_start))
    print("Result difference: %f" % (np.abs(AB-tAB).max(), ))
else:
    from keras.models import Sequential
    from keras.layers import Dense, Dropout, Activation
    from keras.optimizers import SGD
    import numpy as np

    data_dim = 20
    nb_classes = 4

    model = Sequential()

    # Dense(64) is a fully-connected layer with 64 hidden units.
    # in the first layer, you must specify the expected input data shape:
    # here, 20-dimensional vectors.
    model.add(Dense(64, input_dim=data_dim, init='uniform'))
    model.add(Activation('tanh'))
    model.add(Dropout(0.5))
    model.add(Dense(64, init='uniform'))
    model.add(Activation('tanh'))
    model.add(Dropout(0.5))
    model.add(Dense(nb_classes, init='uniform'))
    model.add(Activation('softmax'))

    model.compile(loss='categorical_crossentropy',
                  optimizer='sgd',
                  metrics=["accuracy"])

    # generate dummy training data
    x_train = np.random.random((1000, data_dim))
    y_train = np.random.random((1000, nb_classes))

    # generate dummy test data
    x_test = np.random.random((100, data_dim))
    y_test = np.random.random((100, nb_classes))

    model.fit(x_train, y_train,
              nb_epoch=20,
              batch_size=16)

    score = model.evaluate(x_test, y_test, batch_size=16)


if False:
kernel = conv_peak_shape_kernel(data[:,0], 0, 15)
    cdata=cconv(data[:,1],kernel)
else:
    x = T.ftensor4()
    y = T.ftensor4()
    veclen = x.shape[1]

    conv2d = T.nnet.conv2d

    conv1d_expr = conv2d(x, y, border_mode='full')
    conv1d_expr = theano.sandbox.cuda.dnn.dnn_conv(x,y,border_mode="full")

    conv1d = theano.function([x, y], outputs=conv1d_expr)

    xvals=data[:,1].astype(np.float32)
    fvals=kernel.astype(np.float32)

    cdata = conv1d([[[xvals]]],[[[fvals]]])

print cdata.shape
cdata = np.array(cdata[0,0,0])

print cdata.shape,data.shape
pad=len(kernel)-len(cdata)%len(kernel)
print pad
vpad=np.concatenate([cdata,[0]])
vpad2=np.sum(vpad.reshape(len(kernel),2),axis=1)
print vpad.shape,vpad2.shape
vpad2/=np.amax(vpad2)


'''
