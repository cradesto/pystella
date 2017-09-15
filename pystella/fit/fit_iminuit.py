from iminuit import Minuit


def f(x, y, z):
    return (x - 2) ** 2 + (y - 3) ** 2 + (z - 4) ** 2


m = Minuit(f)

m.migrad()

print(m.values)  # {'x': 2,'y': 3,'z': 4}
print(m.errors)  # {'x': 1,'y': 1,'z': 1}
