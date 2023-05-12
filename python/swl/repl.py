import traceback

import swl.syntax.parser as pr
import swl.syntax.lexer as lr

def repl():
    p = pr.Parser()
    while True:
        try:
            s = input('> ')
            print('tokens:', [t for t in lr.Lexer(s)])
            print('tree:', p.parse(s))
        except EOFError:
            break
        except KeyboardInterrupt:
            break
        except Exception:
            traceback.print_exc()

if __name__ == '__main__':
    repl()

