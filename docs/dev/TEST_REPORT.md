# TEST_REPORT — Relatório de testes e execuções

## Documentação final da Figura 3

Data: 2026-05-06 07:10:45 -03

Objetivo:

```text
Fechar a comparação qualitativa da Figura 3 com a nova política Dirichlet,
decidindo o caso principal e documentando divergências quantitativas.
```

Comandos:

```bash
python3 scripts/plot_cube_plane.py --input data/output/cube_plane_x_5_33_nonuniform.csv --output-dir data/output/figures_nonuniform --all
python3 - <<'PY'  # gera data/output/figure3_region_error_summary.csv
git diff --check
```

Saída do plot da nuvem não uniforme:

```text
input CSV: data/output/cube_plane_x_5_33_nonuniform.csv
points: 10201
grid shape: 101 y samples x 101 z samples
V_num contour: data/output/figures_nonuniform/cube_plane_x_5_33_V_num_contour.png
V_exact contour: data/output/figures_nonuniform/cube_plane_x_5_33_V_exact_contour.png
abs_error contour: data/output/figures_nonuniform/cube_plane_x_5_33_abs_error_contour.png
comparison: data/output/figures_nonuniform/cube_plane_x_5_33_comparison.png
contour levels: 1,2,3,4,5,6,7,8,9,10
V_num isolines: data/output/figures_nonuniform/cube_plane_x_5_33_V_num_isolines.png
V_exact isolines: data/output/figures_nonuniform/cube_plane_x_5_33_V_exact_isolines.png
overlay isolines: data/output/figures_nonuniform/cube_plane_x_5_33_overlay_isolines.png
article style: data/output/figures_nonuniform/cube_plane_x_5_33_article_style.png
metrics: data/output/figures_nonuniform/cube_plane_x_5_33_metrics.txt
metrics CSV: data/output/figures_nonuniform/cube_plane_x_5_33_metrics.csv
plane max abs error: 8.9756567595061476
plane mean abs error: 0.048818265944699984
plane relative error: 0.076503096399384907
interior max abs error: 4.57797229240151
interior mean abs error: 0.042779335597525864
interior relative error: 0.051413379499960249
```

Resumo regional salvo em:

```text
data/output/figure3_region_error_summary.csv
```

```text
regular_refine15,plane_interior,9801,4.5832172421960662,0.033991129037661515,0.045519082240473685
regular_refine15,upper_edge_band_open,200,4.5832172421960662,0.66972900329403773,0.18909845121416047
regular_refine15,central_window,3721,0.068219627386034709,0.018392758048930504,0.0093039163608956994
nonuniform_article_cloud,plane_interior,9801,4.57797229240151,0.042779335597525919,0.051413379499960249
nonuniform_article_cloud,upper_edge_band_open,200,4.57797229240151,0.80546146720387524,0.2075997877939339
nonuniform_article_cloud,central_window,3721,0.090693462740045305,0.019416315028840336,0.0099193017249255763
```

Decisão documentada:

```text
A Figura 3 principal deve usar regular_refine15. A nuvem nonuniform fica como
comparação suplementar: ela tem erro 3D interior próximo, mas erro de plano
maior e GMRES muito mais caro.
```

Conclusão:

```text
PASSOU
```

---

## Política Dirichlet compatível com a série analítica

Data: 2026-05-06 07:03:18 -03

Mudança validada:

```text
As arestas e cantos superiores herdam V = 0 das paredes laterais.
V = V0 é aplicado somente no interior aberto da face z = L.
```

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case sanity
./build/reproduce_cube_sparse --case target
./build/reproduce_cube_sparse --case refine13
./build/reproduce_cube_sparse --case refine15
./build/reproduce_cube_sparse --case plane15
./build/reproduce_cube_sparse --case nonuniform
python3 scripts/plot_cube_plane.py --all
git diff --check
```

Build:

```text
[100%] Built target test_gmres
```

Testes:

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.07 sec
```

Resumo dos casos:

| Caso | nodes | constraints | GMRES | rel_error_global | rel_error_interior | max_abs_error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| sanity 5x5x5 | 125 | 98 | 90 | 3.523634e-01 | 1.894514e-01 | 9.003920e+00 |
| target 11x11x11 | 1331 | 602 | 67 | 4.184668e-02 | 5.134600e-02 | 8.523242e-01 |
| refine13 13x13x13 | 2197 | 866 | 66 | 3.005864e-02 | 4.185387e-02 | 7.438111e-01 |
| refine15 15x15x15 | 3375 | 1178 | 69 | 2.365161e-02 | 2.801970e-02 | 4.642650e-01 |
| nonuniform | 1974 | 762 | 391 | 3.323664e-02 | 2.964210e-02 | 8.344899e-01 |

Critérios de parada:

```text
support_lt_4 = 0 em todos os casos executados
support_lt_8 = 0 em todos os casos executados
mls_failures = 0 em todos os casos executados
gmres_converged = YES em todos os casos executados
```

Exportação `plane15`:

```text
csv:                    data/output/cube_plane_x_5_33_refine15.csv
exported points:        10201
valid points:           10201
MLS failures:           0
max abs error:          8.885380e+00
mean abs error:         3.975798e-02
relative error plane:   7.103024e-02
V_num min:              -3.403056e-02
V_num max:              1.047976e+01
V_exact min:            0.000000e+00
V_exact max:            1.185991e+01
```

Figuras e métricas regeneradas:

```text
data/output/figures/cube_plane_x_5_33_V_num_contour.png
data/output/figures/cube_plane_x_5_33_V_exact_contour.png
data/output/figures/cube_plane_x_5_33_abs_error_contour.png
data/output/figures/cube_plane_x_5_33_comparison.png
data/output/figures/cube_plane_x_5_33_V_num_isolines.png
data/output/figures/cube_plane_x_5_33_V_exact_isolines.png
data/output/figures/cube_plane_x_5_33_overlay_isolines.png
data/output/figures/cube_plane_x_5_33_article_style.png
data/output/figures/cube_plane_x_5_33_metrics.txt
data/output/figures/cube_plane_x_5_33_metrics.csv
```

Conteúdo atualizado de `metrics.txt`:

```text
plane_max_abs_error: 8.8853802327955815
plane_mean_abs_error: 0.039757980669987383
plane_relative_error: 0.071030239356406136
interior_points: 9801
interior_max_abs_error: 4.5832172421960653
interior_mean_abs_error: 0.033991129037661473
interior_relative_error: 0.045519082240473685
contour_levels: 1,2,3,4,5,6,7,8,9,10
article_style_figures_count: 4
```

Conclusão:

```text
PASSOU
```

Observação: os resultados anteriores que explicavam `max_abs_error ≈ V0`
pela prioridade da face superior estão obsoletos para a formulação atual.

---

## Nuvem não uniforme `--case nonuniform`

Data: 2026-05-06

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case nonuniform
```

### Build — nonuniform

```text
[100%] Built target test_gmres
```

Conclusão: build completo sem erros.

### Testes — nonuniform

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.12 sec
```

### Execução do caso nonuniform

```text
=== Nonuniform cloud (base=11 top=13 slices=4 z_frac=0.30) ===

--- Problem ---
  nodes:                  1974
  constraints:            762
  augmented size:         2736

--- MLS diagnostic  (cell-centre grid 15x15x15 = 3375 pts) ---
  invalid (active < 4):   0
  suspect (active < 8):   0
  moment matrix failures: 0
  min active nodes:       8
  max active nodes:       53
  mean active nodes:      24.656888888888886
  max cond estimate:      1.303054e+03

--- Sparse assembly ---
  K nnz  (|v| > 0):       239982
  A_aug COO entries:      269880
  A_aug CSR nnz:          269880
  assembly time:          0.521 s

--- GMRES ---
  restart:                300
  max iter:               20000
  tolerance:              1.0e-09
  iterations:             428
  residual init:          1.676305e+02
  residual final:         1.606441e-07
  rel residual:           9.580e-10
  converged:              YES
  solution time:          0.164 s

--- Errors  (sample grid 11x11x11, valid=1331 interior=729) ---
  max abs error:          1.000000e+01
  rel error (global):     5.529938e-01
  rel error (interior):   4.631267e-02

--- Plane export ---
  csv:                    data/output/cube_plane_x_5_33_nonuniform.csv
  x plane:                5.3300000000000001
  y samples:              101
  z samples:              101
  exported points:        10201
  valid points:           10201
  MLS failures:           0
```

### Critérios de aceitação — nonuniform

```text
support_lt_4:    0   (obrigatório — PASSOU)
support_lt_8:    0   (PASSOU)
mls_failures:    0   (PASSOU)
gmres_converged: YES (PASSOU)
rel error interior: 4,63 % < 5,55 % do 11×11×11 regular (PASSOU)
```

Conclusão:

```text
PASSOU
```

---

## Ajuste de contornos para Figura 3

Data: 2026-05-05 23:01:07 -03

Comandos:

```bash
python3 -m py_compile scripts/plot_cube_plane.py
python3 scripts/plot_cube_plane.py --all
git diff --check
```

Saída do script:

```text
input CSV: data/output/cube_plane_x_5_33_refine15.csv
points: 10201
grid shape: 101 y samples x 101 z samples
V_num contour: data/output/figures/cube_plane_x_5_33_V_num_contour.png
V_exact contour: data/output/figures/cube_plane_x_5_33_V_exact_contour.png
abs_error contour: data/output/figures/cube_plane_x_5_33_abs_error_contour.png
comparison: data/output/figures/cube_plane_x_5_33_comparison.png
contour levels: 1,2,3,4,5,6,7,8,9,10
V_num isolines: data/output/figures/cube_plane_x_5_33_V_num_isolines.png
V_exact isolines: data/output/figures/cube_plane_x_5_33_V_exact_isolines.png
overlay isolines: data/output/figures/cube_plane_x_5_33_overlay_isolines.png
article style: data/output/figures/cube_plane_x_5_33_article_style.png
metrics: data/output/figures/cube_plane_x_5_33_metrics.txt
metrics CSV: data/output/figures/cube_plane_x_5_33_metrics.csv
plane max abs error: 9.9999994786795359
plane mean abs error: 0.042493315964941686
plane relative error: 0.081894083339740326
interior max abs error: 4.3657975314819284
interior mean abs error: 0.03367838888223184
interior relative error: 0.045826314858457826
```

Figuras geradas:

```text
data/output/figures/cube_plane_x_5_33_V_num_contour.png
data/output/figures/cube_plane_x_5_33_V_exact_contour.png
data/output/figures/cube_plane_x_5_33_abs_error_contour.png
data/output/figures/cube_plane_x_5_33_comparison.png
data/output/figures/cube_plane_x_5_33_V_num_isolines.png
data/output/figures/cube_plane_x_5_33_V_exact_isolines.png
data/output/figures/cube_plane_x_5_33_overlay_isolines.png
data/output/figures/cube_plane_x_5_33_article_style.png
```

Conteúdo atualizado de `metrics.txt`:

```text
cube_plane_x_5_33_metrics
input_csv: data/output/cube_plane_x_5_33_refine15.csv
output_dir: data/output/figures
points: 10201
grid_y_count: 101
grid_z_count: 101
x_plane_min: 5.3300000000000001
x_plane_max: 5.3300000000000001
y_min: 0
y_max: 10
z_min: 0
z_max: 10
V_num_min: -0.48030322475030801
V_num_max: 10.03264427824312
V_exact_min: 0
V_exact_max: 11.859907587538387
plane_max_abs_error: 9.9999994786795359
plane_mean_abs_error: 0.042493315964941686
plane_relative_error: 0.081894083339740326
interior_points: 9801
interior_max_abs_error: 4.3657975314819284
interior_mean_abs_error: 0.03367838888223184
interior_relative_error: 0.045826314858457826
contour_levels: 1,2,3,4,5,6,7,8,9,10
article_style_figures_count: 4
article_style_figure: data/output/figures/cube_plane_x_5_33_V_num_isolines.png
article_style_figure: data/output/figures/cube_plane_x_5_33_V_exact_isolines.png
article_style_figure: data/output/figures/cube_plane_x_5_33_overlay_isolines.png
article_style_figure: data/output/figures/cube_plane_x_5_33_article_style.png
```

Conclusão:

```text
PASSOU
```

Observação: build e CTest não foram repetidos nesta rodada porque a alteração
ficou restrita ao script Python de visualização.

## Visualização do plano x = 5.33

Data: 2026-05-05 22:39:02 -03

Comando:

```bash
python3 scripts/plot_cube_plane.py
```

Resultado:

```text
input CSV: data/output/cube_plane_x_5_33_refine15.csv
points: 10201
grid shape: 101 y samples x 101 z samples
V_num contour: data/output/figures/cube_plane_x_5_33_V_num_contour.png
V_exact contour: data/output/figures/cube_plane_x_5_33_V_exact_contour.png
abs_error contour: data/output/figures/cube_plane_x_5_33_abs_error_contour.png
comparison: data/output/figures/cube_plane_x_5_33_comparison.png
metrics: data/output/figures/cube_plane_x_5_33_metrics.txt
metrics CSV: data/output/figures/cube_plane_x_5_33_metrics.csv
plane max abs error: 9.9999994786795359
plane mean abs error: 0.042493315964941686
plane relative error: 0.081894083339740326
interior max abs error: 4.3657975314819284
interior mean abs error: 0.03367838888223184
interior relative error: 0.045826314858457826
```

Arquivos gerados:

```text
data/output/figures/cube_plane_x_5_33_V_num_contour.png
data/output/figures/cube_plane_x_5_33_V_exact_contour.png
data/output/figures/cube_plane_x_5_33_abs_error_contour.png
data/output/figures/cube_plane_x_5_33_comparison.png
data/output/figures/cube_plane_x_5_33_metrics.txt
data/output/figures/cube_plane_x_5_33_metrics.csv
```

Conteúdo do arquivo de métricas:

```text
cube_plane_x_5_33_metrics
input_csv: data/output/cube_plane_x_5_33_refine15.csv
output_dir: data/output/figures
points: 10201
grid_y_count: 101
grid_z_count: 101
x_plane_min: 5.3300000000000001
x_plane_max: 5.3300000000000001
y_min: 0
y_max: 10
z_min: 0
z_max: 10
V_num_min: -0.48030322475030801
V_num_max: 10.03264427824312
V_exact_min: 0
V_exact_max: 11.859907587538387
plane_max_abs_error: 9.9999994786795359
plane_mean_abs_error: 0.042493315964941686
plane_relative_error: 0.081894083339740326
interior_points: 9801
interior_max_abs_error: 4.3657975314819284
interior_mean_abs_error: 0.03367838888223184
interior_relative_error: 0.045826314858457826
```

Conclusão:

```text
PASSOU
```

Observação: build e CTest não foram repetidos nesta rodada porque não houve
alteração em solver, montagem, GMRES ou app C.

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
