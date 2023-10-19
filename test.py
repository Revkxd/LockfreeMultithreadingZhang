# Change the value of these 3 variables to see how Line 7 works
val = 2
thing = 1
zz = 20

ans = 0
ans = 2 if val == 2 else 1 if thing == 1 else 20 if zz == 20 else 0

# Python match statement (similar usage to switch case)
value = 2.0
match value:
    case 1.0: print("value is 1.0")
    case 2.0: print("value is 2.0")
    case _: print("value is neither 1.0 nor 2.0")