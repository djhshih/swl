from swl.transpile._cli import run
from swl.transpile.smk.emit import transpile_dag_file


def main():
    run('Transpile compiled DAG JSON to Snakemake', 'smk', transpile_dag_file, False)
