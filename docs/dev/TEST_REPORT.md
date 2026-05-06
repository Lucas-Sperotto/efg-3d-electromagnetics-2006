# TEST_REPORT — Relatório de testes e execuções

## Análise regional da norma integral da Eq. (16)

Data: 2026-05-06

Objetivo:

```text
Usar a solução GMRES do caso regular refine15, com montagem Gauss ordem 2, e
avaliar a norma integral da Eq. (16) em sub-regiões do cubo usando norm_order=6.
```

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case refine15 --error-region-study
git diff --check
```

Build e testes:

```text
100% tests passed, 0 tests failed out of 29
```

CSV gerado:

- `data/output/error_region_integral_refine15.csv`

Resultados:

| Região | Pontos usados | Eq.16 integral | Numerador | Fração do numerador global | Mean abs integral | Max abs sampled | MLS failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| full_domain | 729000 | 9.180453e-02 | 7.152526e+01 | 100.00 % | 5.067690e-02 | 9.651753e+00 | 0 |
| open_interior | 729000 | 9.180453e-02 | 7.152526e+01 | 100.00 % | 5.067690e-02 | 9.651753e+00 | 0 |
| core_delta_0_5 | 551368 | 2.619258e-02 | 3.453981e+00 | 4.83 % | 2.895232e-02 | 8.187076e-01 | 0 |
| core_delta_1_0 | 373248 | 1.532545e-02 | 6.940169e-01 | 0.97 % | 2.199722e-02 | 3.474524e-01 | 0 |
| central_box | 157464 | 1.027466e-02 | 9.945289e-02 | 0.14 % | 1.789937e-02 | 9.034382e-02 | 0 |
| upper_layer | 72900 | 1.145645e-01 | 7.026686e+01 | 98.24 % | 3.496991e-01 | 9.651753e+00 | 0 |
| upper_edge_bands | 26244 | 2.574390e-01 | 6.947689e+01 | 97.14 % | 8.339680e-01 | 9.651753e+00 | 0 |
| upper_face_interior | 46656 | 1.354568e-02 | 7.899681e-01 | 1.10 % | 7.729788e-02 | 7.978232e-01 | 0 |

Conclusão:

```text
PASSOU. A região que domina o erro é a camada superior z in [9,10], com 98,24 %
do numerador global. Dentro dela, as bandas de aresta superiores respondem por
97,14 % do numerador global; o interior da face superior contribui só 1,10 %.
```

Observação: `open_interior` coincide numericamente com `full_domain` nesta
quadratura porque os pontos de Gauss ficam no interior das células; as faces de
contorno têm medida volumétrica zero.

---

## Estudo da ordem de integração de Gauss-Legendre

Data: 2026-05-06

Objetivo:

```text
Separar a ordem de quadratura usada na montagem de K da ordem usada na norma
integral da Eq. (16), mantendo ordem 2 como padrão, e verificar numericamente
o caso regular refine15.
```

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case refine15 --integration-study norm
./build/reproduce_cube_sparse --case refine15 --integration-study assembly
```

Build e testes:

```text
100% tests passed, 0 tests failed out of 29
```

Arquivos CSV gerados:

- `data/output/integration_order_norm_sensitivity_refine15.csv`
- `data/output/integration_order_solution_sensitivity_refine15.csv`

Parte A — norma Eq. (16) com solução fixa, montagem ordem 2:

| norm_order | Eq.16 integral | Numerador | Denominador | MLS failures |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 6.666831e-02 | 3.668030e+01 | 8.252660e+03 | 0 |
| 2 | 8.217121e-02 | 5.711177e+01 | 8.458361e+03 | 0 |
| 3 | 8.446984e-02 | 6.048276e+01 | 8.476728e+03 | 0 |
| 4 | 9.275266e-02 | 7.302295e+01 | 8.488025e+03 | 0 |
| 5 | 9.186345e-02 | 7.161767e+01 | 8.486619e+03 | 0 |
| 6 | 9.180453e-02 | 7.152526e+01 | 8.486551e+03 | 0 |
| 7 | 9.186332e-02 | 7.161763e+01 | 8.486637e+03 | 0 |
| 8 | 9.185189e-02 | 7.159968e+01 | 8.486622e+03 | 0 |

Parte B — solução variando ordem da montagem, norma Eq. (16) com ordem 6:

| assembly_order | GMRES | Iter | Rel residual | Eq.16 integral | Discrete global | Discrete interior |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | YES | 113 | 8.786874e-10 | 9.842689e-02 | 3.270941e-02 | 4.658729e-02 |
| 2 | YES | 96 | 7.793074e-10 | 9.180453e-02 | 2.365161e-02 | 2.801970e-02 |
| 3 | YES | 96 | 7.345654e-10 | 9.201634e-02 | 2.384583e-02 | 2.846281e-02 |
| 4 | YES | 96 | 8.348048e-10 | 9.205168e-02 | 2.395021e-02 | 2.869954e-02 |
| 5 | YES | 96 | 7.806618e-10 | 9.203966e-02 | 2.388758e-02 | 2.855601e-02 |

Conclusão:

```text
PASSOU. A norma Eq. (16) estabiliza para ordens altas de avaliação; na montagem,
ordem 2 já está no platô observado para refine15. Ordem 1 é visivelmente pior.
```

Observação: o valor histórico `relative_error_global` com Gauss 2x2x2 continua
disponível como padrão; para estudo científico da norma, a coluna com ordem de
referência 6 mostra o valor absoluto mais estável.

---

## Norma relativa integral da Eq. (16)

Data: 2026-05-06

Objetivo:

```text
Trocar a métrica principal relative_error_global da soma discreta em grade
11x11x11 para a integral de domínio da Eq. (16), usando as mesmas células de
integração 15x15x15 e Gauss 2x2x2 da montagem.
```

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case refine15 --solver gmres
./build/reproduce_cube_sparse --case refine15 --solver lapack-dense
./build/reproduce_cube_sparse --case nonuniform_refine --solver gmres
./build/reproduce_cube_sparse --case nonuniform_refine --solver lapack-dense
```

Build e testes:

```text
100% tests passed, 0 tests failed out of 29
```

Resumo dos casos:

| Caso | Solver | Eq.16 domain | Discrete global | Discrete interior | Numerador Eq.16 | Denominador Eq.16 | Quad pts |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| refine15 | gmres | 8.217121e-02 | 2.365161e-02 | 2.801970e-02 | 5.711177e+01 | 8.458361e+03 | 27000 |
| refine15 | lapack-dense | 8.217121e-02 | 2.365161e-02 | 2.801970e-02 | 5.711177e+01 | 8.458361e+03 | 27000 |
| nonuniform_refine | gmres | 8.035461e-02 | 3.284785e-02 | 2.776750e-02 | 5.461448e+01 | 8.458361e+03 | 27000 |
| nonuniform_refine | lapack-dense | 8.035461e-02 | 3.284785e-02 | 2.776750e-02 | 5.461448e+01 | 8.458361e+03 | 27000 |

Conclusão:

```text
PASSOU. relative_error_global agora é a norma integral de domínio da Eq. (16).
```

Observação: a métrica discreta antiga foi preservada em
`relative_error_discrete_global`; `relative_error_interior` continua sendo uma
métrica discreta diagnóstica.

---

## Solver LAPACK denso opcional para validação de GMRES

Data: 2026-05-06

Objetivo:

```text
Validar o sistema aumentado esparso [K G^T; G 0] contra LAPACK dgesv,
sem remover o GMRES atual nem alterar a montagem esparsa.
```

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case refine15 --solver gmres
./build/reproduce_cube_sparse --case refine15 --solver lapack-dense
./build/reproduce_cube_sparse --case nonuniform_refine --solver gmres
./build/reproduce_cube_sparse --case nonuniform_refine --solver lapack-dense
```

Build e testes:

```text
LAPACK found: enabling optional lapack-dense solver
100% tests passed, 0 tests failed out of 29
```

Resumo dos casos:

| Caso | Solver | Iter GMRES | LAPACK info | Residual final | Dif. GMRES/LAPACK | Erro interior 3D | Erro plano |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| refine15 | gmres | 96 | -1 | 1.013100e-07 | n/a | 2.801970e-02 | n/a |
| refine15 | lapack-dense | 96 | 0 | 1.722711e-13 | 8.285314e-10 | 2.801970e-02 | n/a |
| nonuniform_refine | gmres | 858 | -1 | 7.799122e-08 | n/a | 2.776750e-02 | 6.801633e-02 |
| nonuniform_refine | lapack-dense | 858 | 0 | 1.975243e-13 | 3.368672e-08 | 2.776750e-02 | 6.801633e-02 |

Tempos LAPACK:

| Caso | CSR -> dense | dgesv fator/solve | Total solver |
| --- | ---: | ---: | ---: |
| refine15 | 0.108324 s | 1.225903 s | 1.334227 s |
| nonuniform_refine | 0.112976 s | 0.979491 s | 1.092467 s |

Conclusão:

```text
PASSOU. LAPACK e GMRES ficaram próximos nos dois casos executados.
```

Observação: `refine15` não exporta plano nessa CLI; o erro de plano é registrado
quando aplicável, como em `nonuniform_refine`.

---

## Nuvem não uniforme refinada `--case nonuniform_refine`

Data: 2026-05-06 07:18:17 -03

Comandos:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --help
./build/reproduce_cube_sparse --case nonuniform_refine
python3 scripts/plot_cube_plane.py --input data/output/cube_plane_x_5_33_nonuniform_refine.csv --output-dir data/output/figures_nonuniform_refine --all
```

Build e testes:

```text
[100%] Built target test_gmres
100% tests passed, 0 tests failed out of 28
```

Saída do caso:

```text
--- Plane export ---
  csv:                    data/output/cube_plane_x_5_33_nonuniform_refine.csv
  exported points:        10201
  valid points:           10201
  MLS failures:           0
  max abs error:          8.531163e+00
  mean abs error:         6.089040e-02
  relative error plane:   6.801632e-02
  V_num min:              -6.190531e-02
  V_num max:              1.060400e+01
  V_exact min:            0.000000e+00
  V_exact max:            1.185991e+01

=== nonuniform refine cloud (base=13 top=15 slices=4 z_frac=0.30) ===

--- Problem ---
  nodes:                  3089
  constraints:            1082
  augmented size:         4171

--- MLS diagnostic ---
  invalid (active < 4):   0
  suspect (active < 8):   0
  moment matrix failures: 0
  min active nodes:       8
  max active nodes:       54
  mean active nodes:      26.82
  max cond estimate:      1.811e+03

--- Sparse assembly ---
  K nnz  (|v| > 0):       400859
  A_aug CSR nnz:          447167
  assembly time:          0.828 s

--- GMRES ---
  iterations:             832
  residual init:          1.700000e+02
  residual final:         1.684110e-07
  rel residual:           9.907e-10
  converged:              YES
  solution time:          0.716 s

--- Errors ---
  max abs error:          9.407446e-01
  rel error (global):     3.284785e-02
  rel error (interior):   2.776750e-02
```

Métricas do plano refinado:

```text
plane relative error:     6.8016324822001886e-02
interior relative error:  4.4404035267538333e-02
interior mean abs error:  5.6391008910398825e-02
central window rel error: 2.8961925176815718e-02
upper edge open rel:      1.7082540663774520e-01
```

Conclusão:

```text
PASSOU, mas com alerta de custo GMRES.
```

Interpretação: o refinamento melhora levemente o erro 3D interior em relação ao
regular `15x15x15` (`2,776750 %` contra `2,801970 %`) e melhora o erro relativo
interior do plano (`4,4404 %` contra `4,5519 %`), mas piora a janela central e
eleva o GMRES para `832` iterações.

---

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
  max abs error:          9.999999e+00  # Este valor alto no plano é explicado pelo fenômeno de Gibbs da série analítica truncada e pela natureza não interpolante do MLS nas regiões de forte gradiente e descontinuidade de contorno. Consulte `docs/figure3_reproduction.md`, Seção 5.

---

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
