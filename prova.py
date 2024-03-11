import pandas as pd

# Create two DataFrames
df1 = pd.DataFrame({'A': ['A0', 'A1', 'A2'],
                    'B': ['B0', 'B1', 'B2']},
                   index=['K0', 'K1', 'K2'])

df2 = pd.DataFrame({'C': ['C0', 'C2', 'C3'],
                    'D': ['D0', 'D2', 'D3']},
                   index=['K0', 'K2', 'K3'])

# Perform a left join
left_join = df1.merge(df2, left_index=True, right_index=True, how='left')

# Perform an outer join
outer_join = df1.merge(df2, left_index=True, right_index=True, how='outer')

print("Left Join:\n", left_join)
print("\nOuter Join:\n", outer_join)