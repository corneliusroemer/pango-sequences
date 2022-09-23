#%%
import typer 
import pandas as pd
from pango_aliasor.aliasor import Aliasor
import ipdb

def main(
    input_file: str = typer.Argument(..., help="Input file"),
    output_file: str = typer.Argument(..., help="Output file"),
):
    # Load the data
    df = pd.read_csv(input_file, sep="\t")

    field_counter = 10
    for field in ["substitutions", "aaSubstitutions", "deletions", "aaDeletions"]:
        new = []
        new_reversions = []
        for unaliased in df["unaliased"].values:
            parent_unaliased = ".".join(unaliased.split(".")[0:-1])
            parent = df[df["unaliased"] == parent_unaliased]
            this = df[df["unaliased"] == unaliased]
            if unaliased == parent_unaliased:
                new.append("")
                new_reversions.append("")
            else:
                try:
                    parent_nucs = set(parent[field].values[0].split(","))
                except:
                    parent_nucs = set()
                try:
                    this_nucs = set(this[field].values[0].split(","))
                except:
                    this_nucs = set()
            # new.append(",".join(sorted(list(this_nucs - parent_nucs),key=lambda x: int(x[1:-1]))))
            # new_reversions.append(",".join(sorted(list(parent_nucs - this_nucs),key=lambda x: int(x[1:-1]))))
            new.append(",".join(this_nucs - parent_nucs))
            new_reversions.append(",".join(parent_nucs - this_nucs))
        df.insert(field_counter, column=f"new_{field}", value=new)
        df.insert(field_counter, column=f"new_{field}_reversions", value=new_reversions)
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)