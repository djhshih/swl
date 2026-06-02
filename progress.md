
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

