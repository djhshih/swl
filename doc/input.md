
## Workflow with `map` or `map_by`

swl expects all inputs to a record, so the required input for a workflow with
`map` or `map_by` is a record of arrays.

Nextflow and WDL expects the data to be an array of records, so the input data
need to be re-organized for the transpiled workflows for these languages.

