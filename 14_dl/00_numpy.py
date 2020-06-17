import numpy as np

EPOCHS = 300
M = 64
N_I = 1000
N_H = 100
N_O = 10
LEARNING_RATE = 1.0e-06

# create random input and output data
x = np.random.randn(M, N_I)
y = np.random.randn(M, N_O)

# randomly initialize weights
w1 = np.random.randn(N_I, N_H)
w2 = np.random.randn(N_H, N_O)

for t in range(EPOCHS):
    # forward pass
    h = x.dot(w1) # h = x * w1
    h_r = np.maximum(h, 0) # h_r = ReLU(h)
    y_p = h_r.dot(w2) # y_p = h_r * w2

    # compute mean squared error and print loss
    loss = np.square(y_p - y).sum()
    print(t, loss)

    # backward pass: compute gradients of loss with respect to w2
    grad_y_p = 2.0 * (y_p - y)
    grad_w2 = h_r.T.dot(grad_y_p)

    # backward pass: compute gradients of loss with respect to w1
    grad_h_r = grad_y_p.dot(w2.T)
    grad_h = grad_h_r.copy()
    grad_h[h < 0] = 0
    grad_w1 = x.T.dot(grad_h)

    # update weights
    w1 -= LEARNING_RATE * grad_w1
    w2 -= LEARNING_RATE * grad_w2
