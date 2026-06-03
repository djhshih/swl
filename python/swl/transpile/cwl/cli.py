from swl.transpile._cli import run
from swl.transpile.cwl.emit import transpile_dag_file


def main():
    run('Transpile compiled DAG JSON to packed CWL', 'cwl', transpile_dag_file, True)
