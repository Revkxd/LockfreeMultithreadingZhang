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

flatten_list = lambda y:[x for a in y for x in flatten_list(a)] if type(y) is list else [y]
test_arr = [[0.0, 0.0, 10.0], [0.0, 5.0, 0.0]]
print(flatten_list(test_arr))