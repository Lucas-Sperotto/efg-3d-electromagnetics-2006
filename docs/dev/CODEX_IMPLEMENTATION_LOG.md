# CODEX_IMPLEMENTATION_LOG — Registro de implementação

## 2026-05-06 — Nuvem não uniforme `cube_generate_article_cloud`

Objetivo da rodada:

```text
Implementar nuvem não uniforme inspirada na Fig. 2 do artigo (Parreira et al., 2006),
usando grade base regular reforçada com fatias interiores na zona superior do cubo,
e executar o caso nonuniform para comparar com os resultados regulares.
```

Arquivos modificados:

- `src/cube_problem.c` — duas novas funções: `cube_article_cloud_max_nodes` e `cube_generate_article_cloud`
- `include/cube_problem.h` — declarações das duas novas funções
- `apps/reproduce_cube_sparse.c` — enum `CUBE_SPARSE_NONUNIFORM`, função `run_case_nonuniform`, CLI `--case nonuniform`

Arquivos de saída gerados (não versionados):

- `data/output/cube_plane_x_5_33_nonuniform.csv` — plano x = 5.33 da nuvem não uniforme (10 201 pontos)
- `data/output/reproduce_cube_sparse_report.csv` — linha adicionada com resultado nonuniform

Resumo técnico:

**Diagnóstico da tentativa anterior:**

A função `cube_generate_top_refined_nodes` (base_n=7, refine_n=11) gerava edge lines
(t, 0, z) e (t, L, z) que são inteiramente nós de superfície (y=0 e y=L são faces do cubo).
Resultado: 786 restrições em 1308 nós (60 % de superfície) → condicionamento degradado →
erro interior 15,3 %.

**Nova estratégia — `cube_generate_article_cloud`:**

1. Grade regular base_n³ como esqueleto.
2. Face z = L com grade top_n × top_n completa (todos Dirichlet V = V₀, desejado).
3. Para s = 1 .. n_extra_slices: z_s = L − s × (L × z_frac / n_extra_slices).
   Adiciona apenas ix = 1..top_n−2, iy = 1..top_n−2 (interior em xy), sem nós laterais extras.
4. Deduplicação via `append_unique_node`.

`cube_article_cloud_max_nodes` retorna base_n³ + top_n² + n_extra_slices × (top_n−2)² como
cota superior conservadora.

**Parâmetros executados:** base_n=11, top_n=13, n_extra_slices=4, z_frac=0.30.

**Resultados:**

| Métrica | Nuvem não uniforme | 11×11×11 regular | 15×15×15 regular |
| --- | --- | --- | --- |
| Nós totais | 1974 | 1331 | 3375 |
| Restrições | 762 (38,6 %) | 602 (45,2 %) | 1178 (34,9 %) |
| support_lt_4 | 0 | 0 | 0 |
| mls_failures | 0 | 0 | 0 |
| GMRES convergiu | YES | YES | YES |
| GMRES iterações | 428 | — | 70 |
| Erro interior | 4,63 % | 5,55 % | 2,81 % |
| max_abs_error | ~V₀ | ~V₀ | ~V₀ |

A nuvem melhorou o erro interior em relação ao 11×11×11 regular com número semelhante de nós.
O GMRES exigiu mais iterações (428 vs. 70) porque os raios de suporte variam de 8 a 53
(vs. 8 a 27 no regular), tornando o sistema menos uniformemente condicionado.

Comandos executados:

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --case nonuniform
```

Saída principal:

```text
=== Nonuniform cloud (base=11 top=13 slices=4 z_frac=0.30) ===

--- Problem ---
  nodes:                  1974
  constraints:            762
  augmented size:         2736

--- MLS diagnostic ---
  invalid (active < 4):   0
  suspect (active < 8):   0
  moment matrix failures: 0
  min active nodes:       8
  max active nodes:       53
  mean active nodes:      24.66
  max cond estimate:      1303.1

--- Sparse assembly ---
  K nnz  (|v| > 0):       239982
  A_aug CSR nnz:          269880
  assembly time:          0.521 s

--- GMRES ---
  restart:                300
  max iter:               20000
  tolerance:              1.0e-09
  iterations:             428
  residual init:          1.676302e+02
  residual final:         1.606441e-07
  rel residual:           9.58e-10
  converged:              YES
  solution time:          0.164 s

--- Errors (sample grid 11x11x11) ---
  max abs error:          1.000000e+01
  rel error (global):     5.530e-01
  rel error (interior):   4.631e-02
```

Problemas encontrados:

- Nenhum. Build limpo, 28/28 testes passando, GMRES convergido, support_lt_4 = 0.

Próxima etapa recomendada:

```text
Ajustar z_frac ou n_extra_slices para concentrar mais nós na zona de gradiente crítico
(y → 0 ou y → L, z → L), ou implementar precondicionador Jacobi para reduzir iterações GMRES.
```

---

## 2026-05-05 23:01:07 -03 — Iso-linhas para comparação com Figura 3

Objetivo da rodada:

```text
Ajustar apenas scripts/plot_cube_plane.py para gerar figuras mais próximas da
Figura 3 do artigo, usando iso-linhas com níveis fixos de potencial.
```

Arquivos modificados:

- `scripts/plot_cube_plane.py`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `docs/dev/TEST_REPORT.md`
- `TODO.md`

Arquivos gerados em `data/output/figures/`:

- `cube_plane_x_5_33_V_num_contour.png`
- `cube_plane_x_5_33_V_exact_contour.png`
- `cube_plane_x_5_33_abs_error_contour.png`
- `cube_plane_x_5_33_comparison.png`
- `cube_plane_x_5_33_V_num_isolines.png`
- `cube_plane_x_5_33_V_exact_isolines.png`
- `cube_plane_x_5_33_overlay_isolines.png`
- `cube_plane_x_5_33_article_style.png`
- `cube_plane_x_5_33_metrics.txt`
- `cube_plane_x_5_33_metrics.csv`

Resumo técnico:

- Preservadas as figuras preenchidas existentes (`contourf`) e a comparação
  em três painéis.
- Adicionadas iso-linhas com níveis fixos:

  ```text
  1,2,3,4,5,6,7,8,9,10
  ```

- Adicionados modos CLI:
  - `python3 scripts/plot_cube_plane.py`;
  - `python3 scripts/plot_cube_plane.py --article-style`;
  - `python3 scripts/plot_cube_plane.py --all`.
- No modo `--all`, o script gera:
  - iso-linhas de `V_num`;
  - iso-linhas de `V_exact`;
  - sobreposição com `V_num` sólido e `V_exact` tracejado;
  - figura estilo artigo em fundo branco, sem preenchimento.
- O arquivo de métricas registra os níveis de contorno e as figuras estilo
  artigo geradas.
- Nenhuma alteração foi feita em solver, montagem esparsa, GMRES, CSV exportado
  ou `apps/reproduce_cube_sparse.c`.

Comandos executados:

```bash
python3 -m py_compile scripts/plot_cube_plane.py
python3 scripts/plot_cube_plane.py --all
git diff --check
```

Saída principal:

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

Problemas encontrados:

- Nenhum erro de sintaxe no script.
- Nenhuma alteração necessária no pipeline numérico.
- A comparação visual ainda deve ser documentada com a limitação de que a
  malha atual é regular, enquanto o artigo usa nuvem não uniforme.

Próxima etapa recomendada:

```text
Documentar a Figura 3 usando a figura overlay/article_style como comparação
principal, citando os níveis 1..10 V e as limitações de grade regular.
```

---

## 2026-05-05 22:39:02 -03 — Visualização do plano x = 5.33

Objetivo da rodada:

```text
Criar script Python para visualizar o CSV do plano x = 5.33 e gerar figuras
de contorno comparáveis à Figura 3 do artigo.
```

Arquivos modificados:

- `scripts/plot_cube_plane.py`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `docs/dev/TEST_REPORT.md`
- `TODO.md`

Arquivos gerados em `data/output/figures/`:

- `cube_plane_x_5_33_V_num_contour.png`
- `cube_plane_x_5_33_V_exact_contour.png`
- `cube_plane_x_5_33_abs_error_contour.png`
- `cube_plane_x_5_33_comparison.png`
- `cube_plane_x_5_33_metrics.txt`
- `cube_plane_x_5_33_metrics.csv`

Resumo técnico:

- Adicionado `scripts/plot_cube_plane.py`.
- O script usa `pandas`, `numpy` e `matplotlib`, sem `seaborn`.
- O script valida existência do CSV e colunas obrigatórias:
  `x,y,z,V_num,V_exact,abs_error`.
- Os dados são reorganizados em grade regular `y-z` com `y` no eixo horizontal
  e `z` no eixo vertical.
- São geradas três figuras individuais de contorno e uma figura comparativa.
- Na figura comparativa, `V_num` e `V_exact` usam a mesma escala de cor.
- As métricas do plano completo e do interior, com `0 < y < L` e `0 < z < L`,
  são salvas em `cube_plane_x_5_33_metrics.txt`.
- As mesmas métricas também são salvas em CSV para atender à regra geral de
  saída numérica do projeto.
- Nenhuma alteração foi feita em solver, montagem esparsa, GMRES ou
  `apps/reproduce_cube_sparse.c`.

Comandos executados:

```bash
python3 -m py_compile scripts/plot_cube_plane.py
python3 scripts/plot_cube_plane.py
```

Saída principal:

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

Problemas encontrados:

- Nenhuma falha na leitura do CSV.
- Nenhuma lacuna na grade `101x101`.
- O erro máximo do plano permanece dominado pelas bordas, por isso as métricas
  internas foram registradas separadamente.

Próxima etapa recomendada:

```text
Comparar visualmente os contornos gerados com a Figura 3 do artigo e, se
necessário, ajustar apenas o script de visualização para níveis de contorno
iguais aos da publicação.
```

---

## Data

2026-05-05 22:24:42 -03

## Objetivo da rodada

Exportar dados comparáveis ao plano `x = 5.33` do artigo para o caso regular
15x15x15, sem alterar a formulação numérica já validada.

## Arquivos modificados

- `apps/reproduce_cube_sparse.c`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `docs/dev/TEST_REPORT.md`
- `TODO.md`

## Arquivos criados

- Nenhum arquivo versionado.

## Saída numérica CSV

Gerado:

```text
data/output/cube_plane_x_5_33_refine15.csv
```

O arquivo é ignorado por `.gitignore`.

## Resumo técnico das alterações

- Adicionado modo CLI `--case plane15` ao app `reproduce_cube_sparse`.
- O modo `plane15` resolve o caso regular 15x15x15 com:
  - `L = 10`;
  - `V0 = 10`;
  - integração 15x15x15;
  - quadratura Gauss 2x2x2;
  - suporte pelo terceiro vizinho multiplicado por 1.5;
  - Dirichlet via multiplicadores de Lagrange;
  - montagem esparsa;
  - GMRES `restart=300`, `max_iter=20000`.
- Depois da solução, o app exporta o plano `x = 5.33` em grade `101x101`
  no plano `yz`.
- O CSV contém:
  `x,y,z,V_num,V_exact,abs_error`.
- A exportação imprime métricas do plano:
  pontos exportados, falhas MLS, erro máximo, erro médio, erro relativo do
  plano e faixas de `V_num`/`V_exact`.
- Os casos `sanity`, `target`, `refine13`, `refine15` e `all` foram preservados.

## Comandos executados

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --help
./build/reproduce_cube_sparse --case plane15
head -n 5 data/output/cube_plane_x_5_33_refine15.csv
tail -n 5 data/output/cube_plane_x_5_33_refine15.csv
wc -c data/output/cube_plane_x_5_33_refine15.csv
wc -l data/output/cube_plane_x_5_33_refine15.csv
git diff --check
```

## Resultado do build

```text
[100%] Built target test_gmres
```

Build completo sem erros.

## Resultado dos testes

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.12 sec
```

## Resultado da exportação

```text
csv: data/output/cube_plane_x_5_33_refine15.csv
x plane: 5.3300000000000001
y samples: 101
z samples: 101
exported points: 10201
valid points: 10201
MLS failures: 0
max abs error: 9.999999e+00
mean abs error: 4.249332e-02
relative error plane: 8.189408e-02
V_num min: -4.803032e-01
V_num max: 1.003264e+01
V_exact min: 0.000000e+00
V_exact max: 1.185991e+01
```

O sistema 15x15x15 também manteve:

```text
support_lt_4: 0
support_lt_8: 0
mls_failures: 0
gmres_converged: YES
relative_error_interior: 2.812108e-02
```

## CSV

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

## Problemas encontrados

- Nenhuma falha MLS no plano.
- O erro máximo no plano fica próximo de `V0`, como esperado nas regiões de
  borda/topo onde a condição de contorno é contraditória com a série analítica.
- Não foram gerados gráficos nesta rodada.

## Próxima etapa recomendada

Criar script Python de contorno para ler
`data/output/cube_plane_x_5_33_refine15.csv` e gerar mapas de `V_num`,
`V_exact` e `abs_error`.
