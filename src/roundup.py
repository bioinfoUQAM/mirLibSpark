def roundup(x): 
  return x if x % 10 == 0 else x + 10 - x % 10


x = 105
y = roundup(x)
print(x, y)
