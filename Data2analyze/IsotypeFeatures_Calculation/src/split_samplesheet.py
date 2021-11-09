#!/usr/bin/env python3

import click
import pandas as pd
from pathlib import Path


@click.command()
@click.argument('samplesheet', type=click.Path(exists=True), required=True)
@click.option('--output-dir', '-o', type=click.Path(), required=True)
def main(samplesheet, output_dir):
    """ Split migec samplesheet into chunks, such that no chunk shares
    the same fastqs.
    """

    HEADER = [ 'sample', 'master_barcode_seq', 'slave_barcode_seq', 'read1', 'read2' ]
    df = pd.read_csv(samplesheet, dtype=str, header=None, delimiter='\t', names=HEADER)

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for index, (_, frame) in enumerate(df.groupby(['read1'])):
        output = output_dir.joinpath(f'{index}.tsv')
        frame = frame.fillna('')
        frame.to_csv(output, index=False, header=False, sep='\t')


if __name__ == "__main__":
    main()
