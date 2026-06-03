from swl.transpile._cli import run
from swl.transpile.nf.emit import transpile_dag_file


def main():
    run('Transpile compiled DAG JSON to Nextflow DSL2', 'nf', transpile_dag_file, False)
