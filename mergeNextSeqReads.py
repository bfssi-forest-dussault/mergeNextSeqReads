from pathlib import Path
from subprocess import Popen
import click


def concatenate_reads(read_list: [Path], outpath: Path):
    print(f'Concatenating {",".join([x.name for x in read_list])}')
    read_list_text = " ".join([str(read) for read in read_list])
    cmd = f"cat {read_list_text} > {outpath}"
    run_subprocess(cmd)


def run_subprocess(cmd):
    p = Popen(cmd, shell=True)
    p.wait()


def process_sample_dict(sample_dict: dict, outdir: Path):
    for sample_id, reads in sample_dict.items():
        fwd_out = outdir / f'{sample_id}_R1.fastq.gz'
        rev_out = outdir / f'{sample_id}_R2.fastq.gz'
        concatenate_reads(reads['fwd'], fwd_out)
        concatenate_reads(reads['rev'], rev_out)


def group_samples(sample_dir: Path) -> dict:
    """
    Groups samples according sample ID and returns a dict with sample IDs as keys, grouping data across lanes and across
    read direction (R1, R2)

    WARNING: This expects sample file names to follow the BMH naming conventions, e.g.
    BMH-2021-000837_L001_R2.fastq.gz

    """
    samples = list(sample_dir.glob("*.fastq.gz"))
    sample_ids = list(set([x.name.split("_")[0] for x in samples]))
    merge_dict = {}
    for s in sample_ids:
        merge_dict[s] = {'fwd': list(sample_dir.glob(f"{s}*R1.fastq.gz")),
                         'rev': list(sample_dir.glob(f"{s}*R2.fastq.gz"))}
    return merge_dict


def verify_merge_dict(merge_dict: dict):
    for sample_id, reads in merge_dict.items():
        assert len(reads['fwd'] == 4)
        assert len(reads['rev'] == 4)


@click.command(help="Given and input directory containing .fastq.gz files produced by a NextSeq, will merge across all "
                    "4 sample lanes for each sample and deposit the output into a new directory.\n"
                    "WARNING: This software expects the input directory to contain .fastq.gz files that follow the "
                    "standard BMH naming convention, e.g. 'BMH-2021-000837_L001_R2.fastq.gz'")
@click.option('-i', '--indir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the directory on BaseMount for a particular Project. e.g. ../basemount/Projects/[project]. '
                   'Cannot be used at the same time as the --run-dir parameter.')
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Directory to store output merged .fastq.gz files')
def cli(indir: Path, outdir: Path):
    outdir = Path(outdir)
    indir = Path(indir)
    outdir.mkdir(exist_ok=True)
    assert indir.exists()

    merge_dict = group_samples(indir)
    verify_merge_dict(merge_dict)
    process_sample_dict(merge_dict, outdir)


if __name__ == "__main__":
    cli()
