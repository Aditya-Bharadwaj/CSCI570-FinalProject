import random

LETTERS = ['A', 'G', 'T', 'C']
x_string = ''
y_string = ''
for i in range(0, 10000):
    x_index = random.randint(0, 3)
    y_index = random.randint(0, 3)
    x_string += LETTERS[x_index]
    y_string += LETTERS[y_index]

print(x_string)
print()
print(y_string)
