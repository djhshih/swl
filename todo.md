# TODO

- Add stronger semantic tests for concrete batch root typing and table-column propagation.
- Add stronger force/DAG tests for logical mapped table sources and round-trip stability.
- Implement `map_by`
- Implement table update semantics (`tab // rec`, `rec // tab`, `tab // tab`).
- Review CLI/debug output for canonical symbolic DAG and validation entry points.
- Allow binding statements in a block to reference each other, but no recursion
- Continue auditing `spec.md`/`new.md` gaps and convert each into focused tests plus narrow fixes.
