- Do not make code changes that are not necessary, including deleting relevant
  comments and changing the order of function definitions
- Avoid import re-exports
- Unit tests in tests/unit/ should parse strings and not rely on external files
  or disk IO, except where a test explicitly exercises fixture loading from the
  top-level tests/ directory.
- Integration tests are done by ./test.sh using files in tests/
