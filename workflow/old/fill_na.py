import pandas as pd
import sys

# Function to fill missing values in a dataframe with values from the previous column
def fill_missing_with_previous(df):
    for col in df.columns[1:]:
        # Iterate through each column
        for i in range(0, len(df)):
            if df.loc[i, col] in {"d__", "p__", "c__", "o__", "f__", "g__", "s__"}:
                # If the value is missing, fill it with the previous column's value
                string = df.loc[i, df.columns.get_loc(col) - 1].split("__")[1]
                df.loc[i, col] = df.loc[i, col] + string
    return df

def main():
    # Read input dataframe from stdin
    try:
        df = pd.read_csv(sys.stdin, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        sys.exit("Error: Empty input dataframe")

    # Fill missing values in the dataframe with values from the previous column
    filled_df = fill_missing_with_previous(df)

    # Output the modified dataframe to stdout
    filled_df.to_csv(sys.stdout, sep='\t', index=False, header=None)

if __name__ == "__main__":
    main()
