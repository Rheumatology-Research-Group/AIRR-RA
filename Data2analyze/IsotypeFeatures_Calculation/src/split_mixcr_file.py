#!/usr/bin/env python3

import csv
import pandas as pd
import click
from pathlib import Path


@click.command()
@click.argument('file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', type=click.Path(), required=True)
def main(file, output_dir):
    """ Split mixcr exported clones output file into 7 chain files """

    df = pd.read_csv(file, delimiter='\t', low_memory=False, quoting=csv.QUOTE_ALL)

    name = Path(file).stem
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for chain in ['TRA', 'TRB', 'TRD', 'TRG', 'IGH', 'IGK', 'IGL']:
        sub = df[df['allJHitsWithScore'].str.match(f'{chain}J')]
        output = output_dir.joinpath(f'{name}__{chain}.tsv')
        sub.to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    main()
