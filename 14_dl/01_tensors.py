import torch

EPOCHS = 300
M = 64
N_I = 1000
N_H = 100
N_O = 10
LEARNING_RATE = 1.0e-06

# create random input and output data
x = torch.randn(M, N_I)
y = torch.randn(M, N_O)

# randomly initialize weights
w1 = torch.randn(N_I, N_H)
w2 = torch.randn(N_H, N_O)

for t in range(EPOCHS):
    # forward pass: compute predicted y
    h = x.mm(w1)
    h_r = h.clamp(min=0)
    y_p = h_r.mm(w2)

    # compute and print loss
    loss = (y_p - y).pow(2).sum().item()
    print(t, loss)

    # backward pass: compute gradients of loss with respect to w2
    grad_y_p = 2.0 * (y_p - y)
    grad_w2 = h_r.t().mm(grad_y_p)

    # backward pass: compute gradients of loss with respect to w1
    grad_h_r = grad_y_p.mm(w2.t())
    grad_h = grad_h_r.clone()
    grad_h[h < 0] = 0
    grad_w1 = x.t().mm(grad_h)

    # update weights
    w1 -= LEARNING_RATE * grad_w1
    w2 -= LEARNING_RATE * grad_w2
