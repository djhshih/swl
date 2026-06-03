from swl.transpile._cli import run
from swl.transpile.wdl.emit import transpile_dag_file


def main():
    run('Transpile compiled DAG JSON to WDL 1.1', 'wdl', transpile_dag_file, False)
