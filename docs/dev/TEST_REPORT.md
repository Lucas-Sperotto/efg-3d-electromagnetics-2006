# TEST_REPORT — Relatório de testes e execuções

## Data

2026-05-05 22:24:42 -03

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
Total Test time (real) = 0.12 sec
```

## CLI

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
  ./build/reproduce_cube_sparse --case plane15
  ./build/reproduce_cube_sparse --case all
```

## Exportação do plano

Comando:

```bash
./build/reproduce_cube_sparse --case plane15
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
plane = x = 5.33
y_samples = 101
z_samples = 101
```

Saída:

```text
--- Plane export ---
  csv:                    data/output/cube_plane_x_5_33_refine15.csv
  x plane:                5.3300000000000001
  y samples:              101
  z samples:              101
  exported points:        10201
  valid points:           10201
  MLS failures:           0
  max abs error:          9.999999e+00
  mean abs error:         4.249332e-02
  relative error plane:   8.189408e-02
  V_num min:              -4.803032e-01
  V_num max:              1.003264e+01
  V_exact min:            0.000000e+00
  V_exact max:            1.185991e+01

=== 15x15x15 plane15 (export x=5.33 plane) ===

--- Problem ---
  nodes:                  3375
  constraints:            1178
  augmented size:         4553
  node grid:              15x15x15
  integration cells:      15x15x15

--- MLS diagnostic  (cell-centre grid 15x15x15 = 3375 pts) ---
  invalid (active < 4):   0
  suspect (active < 8):   0
  moment matrix failures: 0
  min active nodes:       8
  max active nodes:       27
  mean active nodes:      23.56
  max cond estimate:      8.205e+02
  mean cond estimate:     5.042e+02

--- Sparse assembly ---
  K nnz  (|v| > 0):       328509
  A_aug COO entries:       368199
  A_aug CSR nnz:           368199
  assembly time:           0.642 s

--- GMRES ---
  restart:                300
  max iter:               20000
  tolerance:              1.0e-09
  iterations:             70
  residual init:          1.500000e+02
  residual final:         1.354699e-07
  rel residual:           9.031e-10
  converged:              YES
  solution time:          0.025 s

--- Errors  (sample grid 11x11x11, valid=1331 interior=729) ---
  max abs error:          1.000002e+01
  rel error (global):     5.530370e-01
  rel error (interior):   2.812108e-02
```

## CSV gerado

Arquivo:

```text
data/output/cube_plane_x_5_33_refine15.csv
```

Tamanho:

```text
1119344 bytes
10202 linhas incluindo cabeçalho
```

Primeiras 5 linhas:

```text
x,y,z,V_num,V_exact,abs_error
5.3300000000000001,0,0,-4.6981770252577896e-08,0,4.6981770252577896e-08
5.3300000000000001,0,0.10000000000000001,0.00011260152361410542,0,0.00011260152361410542
5.3300000000000001,0,0.20000000000000001,0.00014497959621752784,0,0.00014497959621752784
5.3300000000000001,0,0.29999999999999999,0.00010216755348659567,0,0.00010216755348659567
```

Últimas 5 linhas:

```text
5.3300000000000001,10,9.5999999999999996,3.990924108594033,6.8005564903065482e-15,3.9909241085940264
5.3300000000000001,10,9.6999999999999993,5.5101853986590887,9.1101140713595324e-15,5.5101853986590799
5.3300000000000001,10,9.8000000000000007,7.0240861692584167,1.3331103951700499e-14,7.0240861692584033
5.3300000000000001,10,9.9000000000000004,8.5210718264413963,2.2359023490980203e-14,8.5210718264413732
5.3300000000000001,10,10,9.9999994786795821,4.5288305635599732e-14,9.9999994786795376
```

## Critérios de parada

```text
MLS failures no plano: 0
arquivo CSV: escrito com sucesso
gmres_converged: YES
support_lt_4: 0
support_lt_8: 0
```

Conclusão:

```text
PASSOU
```

## Próximo passo recomendado

```text
Criar script Python de contorno para V_num, V_exact e abs_error.
```
