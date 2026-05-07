## Open semantic questions

- Interpolation values are currently preserved as syntax nodes (`Word`, `Var`, `Expr`) in task defaults. This seems correct, but the exact resolution phase still needs to be defined clearly.

  * The value of Word will be substitued verbatim at compile time, and we can do
    compile time checks to ensure that syntax is valid. Var will be substituted with
    the value it holds at runtime. Expr will be evaluated at runtime and
    substituted with the resulting value.

- It is not yet fully settled whether task output params must always have defaults or whether they may be derived only from body/runtime behavior.

   * Task output params must have defaults. The defaults can have glob *
     pattern.

- `bash.py` exists as an optional analyzer, but current semantics should probably treat annotation metadata as authoritative and body text as opaque.

   * We will use bash.py to ensure that the bash syntax after string
     interpolation is still valid as a pre-runtime and runtime check.

- Tighten semantics around non-record application arguments 
  
   * Applying task to a scalar argument is legal. The scalar value will be lifted
     to a record with a single field that has the same name as the first input
     of the task.

- Clarify partial application semantics 
   
   * Partial application creates a function enclosed with the provided inputs
     (closure)

- Separate “issues” from “hard errors”

   * All issues are errors unless specified otherwise.

- Decide whether workflow signatures should preserve types

   * A workflow must evaluate to a function. The workflow output is the output
     of this function. This function may be an explicit lambda, in which case
     the outputs are explicit defined on the final line of body fo the lambda.
     This function can be a named task, so the workflow output is the same as
     as the task. This function can be the result of a chain, semantic analysis
     should perform an union of the output variables of each task in the chain
     from left to right.
