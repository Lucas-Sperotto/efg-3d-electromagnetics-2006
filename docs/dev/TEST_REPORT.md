# TEST_REPORT — Relatório de testes e execuções

## Data

Preencher ao executar.

## Commit analisado

```bash
git log --oneline -1
````

Resultado:

```text
Preencher.
```

## Estado do Git

```bash
git status --short
```

Resultado:

```text
Preencher.
```

## Build

Comando:

```bash
cmake --build build
```

Resultado:

```text
Preencher.
```

## Testes

Comando:

```bash
ctest --test-dir build --output-on-failure
```

Resultado:

```text
Preencher.
```

## Testes individuais relevantes

* [ ] test_mls_diagnostic
* [ ] test_sparse_matrix
* [ ] test_sparse_global_stiffness
* [ ] test_lagrange_augmented_sparse
* [ ] test_gmres

Resultado:

```text
Preencher.
```

## Execução do caso pequeno

Comando:

```bash
./build/reproduce_cube_sparse
```

Configuração:

```text
nodes:
cells:
quadrature:
solver:
```

Resultado:

```text
Preencher.
```

## Execução do caso 11x11x11

Configuração:

```text
L = 10
V0 = 10
nodes = 11x11x11
cells = 15x15x15
quadrature = 2x2x2
boundary = Dirichlet via Lagrange
solver = GMRES
```

Resultado:

```text
nodes:
constraints:
augmented_size:
K_nnz:
A_aug_nnz:
support_min:
support_mean:
support_max:
support_lt_4:
support_lt_8:
mls_failures:
gmres_converged:
gmres_iterations:
residual_initial:
residual_final:
relative_error_global:
relative_error_internal:
max_abs_error:
assembly_time_s:
solve_time_s:
```

## Interpretação preliminar

### Critérios de bloqueio

```text
support_lt_4:
mls_failures:
gmres_converged:
residual_final:
```

Conclusão:

```text
PASSOU / NÃO PASSOU / INCONCLUSIVO
```

## Próximo passo recomendado

```text
Preencher.
```
