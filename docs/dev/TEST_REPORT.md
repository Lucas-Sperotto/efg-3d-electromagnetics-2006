# TEST_REPORT — Relatório de testes e execuções

## Data

2026-05-05 21:55:38 -03

## Commit analisado

Comando:

```bash
git log --oneline -1
```

Resultado:

```text
fec2d67 audit: complete CLAUDE_AUDIT.md + fix diagnostic cell-centre placement
```

## Estado do Git

Comando:

```bash
git status --short
```

Resultado antes de preencher estes relatórios:

```text
 M apps/reproduce_cube_sparse.c
```

Resultado esperado após esta rodada:

```text
 M apps/reproduce_cube_sparse.c
 M docs/dev/CODEX_IMPLEMENTATION_LOG.md
 M docs/dev/TEST_REPORT.md
```

## Build

Comando:

```bash
cmake --build build
```

Resultado:

```text
[100%] Built target test_gmres
```

Conclusão: build completo sem erros.

## Testes

Comando esperado:

```bash
ctest --test-dir build --output-on-failure
```

Resultado:

```text
Traceback (most recent call last):
  File "/home/sperotto/.local/bin/ctest", line 5, in <module>
    from cmake import ctest
ModuleNotFoundError: No module named 'cmake'
```

Diagnóstico:

```text
command -v ctest
/home/sperotto/.local/bin/ctest
```

Comando efetivo usado para validar a suíte:

```bash
/usr/bin/ctest --test-dir build --output-on-failure
```

Resultado:

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.14 sec
```

## Testes individuais relevantes

```text
test_mls_diagnostic: PASS
test_sparse_matrix: PASS
test_sparse_global_stiffness: PASS
test_lagrange_augmented_sparse: PASS
test_gmres: PASS
```

## Ajuda da CLI

Comando:

```bash
./build/reproduce_cube_sparse --help
```

Resultado:

```text
Usage:
  ./build/reproduce_cube_sparse
  ./build/reproduce_cube_sparse --case sanity
  ./build/reproduce_cube_sparse --case target
  ./build/reproduce_cube_sparse --case all

Cases:
  sanity  regular 5x5x5 nodes, 5x5x5 integration cells, dense comparison
  target  regular 11x11x11 nodes, 15x15x15 integration cells, GMRES
  all     run sanity followed by target (default)
```

## Execução do caso pequeno

Comando:

```bash
./build/reproduce_cube_sparse --case sanity
```

Configuração:

```text
L = 10
V0 = 10
nodes = 5x5x5
cells = 5x5x5
quadrature = 2x2x2
boundary = Dirichlet via Lagrange
solver = GMRES, com comparação densa
```

Resultado:

```text
nodes: 125
constraints: 98
augmented_size: 223
K_nnz: 6859
A_aug_nnz: 9715
support_min: 8
support_mean: 17.576000
support_max: 27
support_lt_4: 0
support_lt_8: 0
mls_failures: 0
gmres_converged: YES
gmres_iterations: 89
residual_initial: 5.000000e+01
residual_final: 4.502732e-08
relative_error_global: 6.556602e-01
relative_error_internal: 2.178039e-01
max_abs_error: 1.000250e+01
assembly_time_s: 0.009975
solve_time_s: 0.000602
GMRES vs dense diff: 1.738e-09
```

## Execução do caso 11x11x11

Comando:

```bash
./build/reproduce_cube_sparse --case target
```

Configuração:

```text
L = 10
V0 = 10
nodes = 11x11x11
cells = 15x15x15
quadrature = 2x2x2
support = terceiro vizinho multiplicado por 1.5
boundary = Dirichlet via Lagrange
assembly = sparse
solver = GMRES
```

Resultado:

```text
nodes: 1331
constraints: 602
augmented_size: 1933
K_nnz: 117649
A_aug_nnz: 137383
support_min: 8
support_mean: 23.557630
support_max: 27
support_lt_4: 0
support_lt_8: 0
mls_failures: 0
gmres_converged: YES
gmres_iterations: 67
residual_initial: 1.100000e+02
residual_final: 8.698756e-08
relative_error_global: 5.533057e-01
relative_error_internal: 5.552359e-02
max_abs_error: 1.000000e+01
assembly_time_s: 0.412413
solve_time_s: 0.008328
```

## Execução padrão

Comando:

```bash
./build/reproduce_cube_sparse
```

Resultado:

```text
Executou os casos sanity e target em sequência.
Ambos retornaram código 0.
Gerou data/output/reproduce_cube_sparse_report.csv.
```

CSV:

```text
label,nodes,constraints,augmented_size,K_nnz,A_aug_nnz,support_min,support_mean,support_max,support_lt_4,support_lt_8,mls_failures,gmres_converged,gmres_iterations,residual_initial,residual_final,relative_error_global,relative_error_internal,max_abs_error,assembly_time_s,solve_time_s
"5x5x5 sanity check (dense comparison enabled)",125,98,223,6859,9715,8,17.575999999999997,27,0,0,0,YES,89,50,4.502731637452715e-08,0.65566020220249721,0.2178038587457618,10.002496316879469,0.0099749999999999995,0.000602
"11x11x11 target (15x15x15 integration cells, GMRES only)",1331,602,1933,117649,137383,8,23.557629629629595,27,0,0,0,YES,67,110,8.698755847818432e-08,0.55330567475497328,0.055523591643745318,10.000000001957561,0.41241299999999997,0.0083280000000000003
```

## Interpretação preliminar

### Critérios de bloqueio

```text
support_lt_4: 0 nos dois casos
mls_failures: 0 nos dois casos
gmres_converged: YES nos dois casos
residual_final:
  sanity: 4.502732e-08
  target: 8.698756e-08
```

Conclusão:

```text
PASSOU
```

O pipeline esparso com Dirichlet via multiplicadores de Lagrange e GMRES está
operacional para o alvo inicial 11x11x11.

## Problemas encontrados

```text
ctest sem caminho usa /home/sperotto/.local/bin/ctest e falha por falta do
módulo Python cmake. /usr/bin/ctest funciona e validou 28/28 testes.
```

## Próximo passo recomendado

```text
Pedir ao Gemini uma revisão focada de apps/reproduce_cube_sparse.c, em especial:
relatório obrigatório, critérios de parada, interpretação de support_lt_8 e
decisão sobre avançar para 13x13x13/15x15x15 antes de exportar o plano x = 5.33.
```
