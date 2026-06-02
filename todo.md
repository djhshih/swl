# TODO

- Add stronger semantic tests for concrete batch root typing and table-column propagation.
- Add stronger force/DAG tests for logical mapped table sources and round-trip stability.
- Implement `map_by`
- Implement table update semantics (`tab // rec`, `rec // tab`, `tab // tab`).
- Review CLI/debug output for canonical symbolic DAG and validation entry points.
- Allow binding statements in a block to reference each other, but no recursion
- Continue auditing `spec.md`/`new.md` gaps and convert each into focused tests plus narrow fixes.

### Partial or in progress

| Feature | Current status |
|---------|----------------|
| Input type/desc materialization | Workflow input metadata is refined from inferred interface information and step specs, but coverage may depend on available inferred schema information. |
| Merge elimination | The specification requires elimination before final DAG emission; implementation work may still be required where transient merges survive too long. |
| Record saturation | Direct-call record saturation is required by the specification; remaining non-saturating cases require explicit record serialization with full field structure. |
| `map.scatter` / `map.broadcast` population | The specification requires these fields in final mapped steps; implementation work may still be needed to guarantee they are always present. |
| Output interface materialization | Final outputs are normalized to explicit top-level outputs; richer output metadata may still depend on available inferred type information. |
| Nested field projection support | Basic field projection is represented in bindings; more complex nested forms may require additional normalization or target-specific handling. |

### Not yet guaranteed everywhere

| Feature | Required contract |
|---------|-------------------|
| Final-DAG merge freedom | Final emitted DAGs must not contain `merge` bindings. |
| Full optionality propagation | Final emitted interface and parameter specs must preserve optionality explicitly. |
| Universal mapped-port classification | Every mapped-step input must appear in exactly one of `map.scatter` or `map.broadcast`. |
| Portable handling of remaining record values | Any `record` that remains in the final DAG must have explicit structure and must not stand in for an unflattened direct step-call argument. |

