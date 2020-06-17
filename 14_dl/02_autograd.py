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
w1 = torch.randn(N_I, N_H, requires_grad=True)
w2 = torch.randn(N_H, N_O, requires_grad=True)

for t in range(EPOCHS):
    # forward pass: compute predicted y
    h = x.mm(w1)
    h_r = h.clamp(min=0)
    y_p = h_r.mm(w2)

    # compute and print loss
    loss = (y_p - y).pow(2).sum()
    print(t, loss.item())

    # backward pass
    loss.backward()

    with torch.no_grad():
        # update weights
        w1 -= LEARNING_RATE * w1.grad
        w2 -= LEARNING_RATE * w2.grad

        # initialize weights
        w1.grad.zero_()
        w2.grad.zero_()
