#%%
import typer 
import pandas as pd
from pango_aliasor.aliasor import Aliasor

def main(
    input_file: str = typer.Argument(..., help="Input file"),
    output_file: str = typer.Argument(..., help="Output file"),
):
    # Initalize aliasor
    aliasor = Aliasor()
    # Load the data
    df = pd.read_csv(input_file, sep="\t")
    unaliased = df["seqName"].apply(aliasor.uncompress)
    df.insert(1, column="unaliased", value=unaliased)
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)