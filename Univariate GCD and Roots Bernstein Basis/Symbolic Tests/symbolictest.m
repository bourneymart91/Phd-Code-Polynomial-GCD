function [] = symbolictest

% initialise x to be symbolic
x = sym('x')
y = sym('y')


a = (1-y)^2 .* x

b = (1-y)*y*(1-x)

c = y*(1-x)^2

d = (1-y)^2 * x

e = (1-y) * y * x

f = y^2 * (1-x)

a + b + c + d + e + f

