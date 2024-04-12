import numpy as np

np.random.seed(33)

# Between -1 and 1 (normal)
def recurs_test(low, high):
    norm_correlation_values = list(np.random.uniform(low = low, high=high, size=100))
    # print(norm_correlation_values)
    norm_naked_sum = 0
    for r in norm_correlation_values:
        w = 2 if r > 0 else 1
        norm_naked_sum += w * (r**2)
    output = f'for low {low} & high {high}: \n result = {norm_naked_sum} \n'
    return output



print(recurs_test(-.25,-.2))
print(recurs_test(0.1,0.2))




# # Between -1 and 0.2 (negative skew). Should give better result than the above 
# neg_correlation_values = list(np.random.uniform(low = -1.0, high=0.2, size=10))
# print(neg_correlation_values)


# neg_naked_sum = 0
# for g in neg_correlation_values:
#     for j in neg_correlation_values:
#         neg_naked_sum =+ g+j


# print(neg_naked_sum)

# # Between -0.2 and 1 (positive skew). Should give worse result than the above and far bett
# pos_correlation_values = list(np.random.uniform(low = -1.0, high=0.2, size=10))
# print(pos_correlation_values)


# pos_naked_sum = 0
# for g in pos_correlation_values:
#     for j in pos_correlation_values:
#         pos_naked_sum =+ g+j


# print(pos_naked_sum)

# # Between -0.5 and 0.5 (tighter). 
# tight_correlation_values = list(np.random.uniform(low = -1.0, high=0.2, size=10))
# print(tight_correlation_values)


# tight_naked_sum = 0
# for g in tight_correlation_values:
#     for j in tight_correlation_values:
#         tight_naked_sum =+ g+j


# print(tight_naked_sum)