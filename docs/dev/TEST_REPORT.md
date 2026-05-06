# TEST_REPORT — Relatório de testes e execuções

## Data

2026-05-05 22:19:00 -03

## Commit analisado

```bash
git log --oneline -1
```

Resultado:

```text
09931b7 feat: validate sparse cube pipeline with GMRES
```

## Estado do Git

```bash
git status --short
```

Resultado antes de atualizar estes relatórios:

```text
 M apps/reproduce_cube_sparse.c
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

Comando:

```bash
/usr/bin/ctest --test-dir build --output-on-failure
```

Resultado:

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.14 sec
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
  ./build/reproduce_cube_sparse --case refine13
  ./build/reproduce_cube_sparse --case refine15
  ./build/reproduce_cube_sparse --case all

Cases:
  sanity  regular 5x5x5 nodes, 5x5x5 integration cells, dense comparison
  target  regular 11x11x11 nodes, 15x15x15 integration cells, GMRES
  refine13 regular 13x13x13 nodes, 15x15x15 integration cells, GMRES
  refine15 regular 15x15x15 nodes, 15x15x15 integration cells, GMRES
  all     run sanity, target, refine13, and refine15 (default)
```

## Execução refine13

Comando:

```bash
./build/reproduce_cube_sparse --case refine13
```

Configuração:

```text
L = 10
V0 = 10
nodes = 13x13x13
cells = 15x15x15
quadrature = 2x2x2
support = terceiro vizinho multiplicado por 1.5
boundary = Dirichlet via Lagrange
assembly = sparse
solver = GMRES restart=300 max_iter=20000
```

Resultado:

```text
nodes: 2197
constraints: 866
augmented_size: 3063
K_nnz: 205379
A_aug_nnz: 234169
support_min: 8
support_mean: 23.557630
support_max: 27
support_lt_4: 0
support_lt_8: 0
mls_failures: 0
max_cond_estimate: 9.004260e+02
gmres_converged: YES
gmres_iterations: 66
residual_initial: 1.300000e+02
residual_final: 9.294019e-08
relative_error_global: 5.535730e-01
relative_error_interior: 4.566296e-02
max_abs_error: 1.000013e+01
assembly_time_s: 0.528650
solve_time_s: 0.016021
```

## Execução refine15

Comando:

```bash
./build/reproduce_cube_sparse --case refine15
```

Configuração:

```text
L = 10
V0 = 10
nodes = 15x15x15
cells = 15x15x15
quadrature = 2x2x2
support = terceiro vizinho multiplicado por 1.5
boundary = Dirichlet via Lagrange
assembly = sparse
solver = GMRES restart=300 max_iter=20000
```

Resultado:

```text
nodes: 3375
constraints: 1178
augmented_size: 4553
K_nnz: 328509
A_aug_nnz: 368199
support_min: 8
support_mean: 23.557630
support_max: 27
support_lt_4: 0
support_lt_8: 0
mls_failures: 0
max_cond_estimate: 8.204838e+02
gmres_converged: YES
gmres_iterations: 70
residual_initial: 1.500000e+02
residual_final: 1.354699e-07
relative_error_global: 5.530370e-01
relative_error_interior: 2.812108e-02
max_abs_error: 1.000002e+01
assembly_time_s: 0.729538
solve_time_s: 0.029133
```

## Execução padrão

Comando:

```bash
./build/reproduce_cube_sparse
```

Resultado: executou `sanity`, `target`, `refine13` e `refine15` em sequência,
retornando código 0.

CSV gerado:

```text
data/output/reproduce_cube_sparse_report.csv
```

Conteúdo:

```text
label,nodes,constraints,augmented_size,K_nnz,A_aug_nnz,support_min,support_mean,support_max,support_lt_4,support_lt_8,mls_failures,max_cond_estimate,gmres_converged,gmres_iterations,residual_initial,residual_final,relative_error_global,relative_error_interior,max_abs_error,assembly_time_s,solve_time_s
"5x5x5 sanity check (dense comparison enabled)",125,98,223,6859,9715,8,17.575999999999997,27,0,0,0,61.413003139977342,YES,89,50,4.502731637452715e-08,0.65566020220249721,0.2178038587457618,10.002496316879469,0.011518,0.00077200000000000001
"11x11x11 target (15x15x15 integration cells, GMRES only)",1331,602,1933,117649,137383,8,23.557629629629595,27,0,0,0,676.05833333380599,YES,67,110,8.698755847818432e-08,0.55330567475497328,0.055523591643745318,10.000000001957561,0.45323000000000002,0.0098860000000000007
"13x13x13 refine13 (15x15x15 integration cells, GMRES only)",2197,866,3063,205379,234169,8,23.557629629629595,27,0,0,0,900.42600000087532,YES,66,130,9.2940185865611113e-08,0.55357298992722137,0.045662955683606739,10.000126758612106,0.49610900000000002,0.013061
"15x15x15 refine15 (15x15x15 integration cells, GMRES only)",3375,1178,4553,328509,368199,8,23.557629629629595,27,0,0,0,820.4837762407044,YES,70,150,1.3546990018799636e-07,0.55303702122071197,0.028121081303306723,10.000018650316381,0.67268099999999997,0.026138000000000002
```

## Comparação 11/13/15

```text
11x11x11: nodes=1331, constraints=602,  A_aug_nnz=137383, interior_error=5.552359e-02, gmres_iterations=67
13x13x13: nodes=2197, constraints=866,  A_aug_nnz=234169, interior_error=4.566296e-02, gmres_iterations=66
15x15x15: nodes=3375, constraints=1178, A_aug_nnz=368199, interior_error=2.812108e-02, gmres_iterations=70
```

O erro interior caiu monotonicamente entre 11, 13 e 15. Não há suspeita de
regressão por erro crescente nesta rodada.

## Critérios de bloqueio

```text
support_lt_4: 0 em 13x13x13 e 15x15x15
support_lt_8: 0 em 13x13x13 e 15x15x15
mls_failures: 0 em 13x13x13 e 15x15x15
gmres_converged: YES em 13x13x13 e 15x15x15
```

Conclusão:

```text
PASSOU
```

## Próximo passo recomendado

```text
Exportar plano x = 5.33 para o caso 15x15x15, com CSV y,z,V_num,V_exact,abs_error.
```
