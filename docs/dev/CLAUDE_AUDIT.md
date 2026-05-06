# CLAUDE_AUDIT — Auditoria matemática e arquitetural

## Data

2026-05-05

## Commit analisado

```
aacbcd0 app: add sparse EFG cube runner with GMRES and MLS diagnostic
```

## Objetivo da auditoria

Avaliar se o projeto está pronto para executar com segurança casos regulares maiores do cubo eletrostático 3D com:

```
L = 10, V0 = 10
grade regular de nós (11x11x11, 13x13x13, 15x15x15)
15x15x15 células de integração
Gauss ordem 2 (2x2x2 por célula)
Dirichlet via multiplicadores de Lagrange
montagem esparsa
solução com GMRES
```

## Arquivos analisados

- [x] `TODO.md`
- [x] `README.md`
- [x] `CMakeLists.txt`
- [x] `include/` (todos os headers)
- [x] `src/mls_diagnostic.c`
- [x] `src/sparse_matrix.c`
- [x] `src/global_stiffness.c`
- [x] `src/lagrange_system.c`
- [x] `src/gmres.c`
- [x] `src/mls.c`
- [x] `src/stiffness.c`
- [x] `src/dirichlet.c`
- [x] `src/analytical.c`
- [x] `src/cube_problem.c`
- [x] `src/support.c`
- [x] `apps/reproduce_cube_sparse.c`
- [x] `tests/test_mls_diagnostic.c`
- [x] `tests/test_sparse_global_stiffness.c`
- [x] `tests/test_gmres.c`
- [x] `docs/02_formulacao_matematica.md`
- [x] `docs/03_resultados_numericos.md`
- [x] `docs/05_mapa_equacoes_codigo.md`
- [x] `docs/06_condicoes_de_contorno_e_lagrange.md`
- [x] `docs/04_parametros_numericos.md`

---

## Estado verificado: 28/28 testes passando

```
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.08 sec
```

---

## 1. Diagnóstico MLS — VERIFICADO CORRETO

### 1.1. Cardinalidade mínima

O arquivo `include/mls_diagnostic.h` define:

```c
#define MLS_MIN_SUPPORT_NODES    4   /* base linear 3D: p = [1, x, y, z] */
#define MLS_SUSPECT_SUPPORT_NODES 8
```

`mls_diagnostic_at_point` classifica corretamente:

- `active_nodes < 4` → `MLS_POINT_INVALID`
- `4 <= active_nodes < 8` → `MLS_POINT_SUSPECT`
- `active_nodes >= 8` → `MLS_POINT_VALID`

### 1.2. Construção da matriz de momento A(x)

`build_moment_matrix` em `src/mls_diagnostic.c`:

```c
A(x) = Σ_I w(x - x_I) p(x_I) p^T(x_I)
```

Usa `nodes[i].x` (não o ponto de avaliação) para a base — correto conforme a formulação MLS.

### 1.3. Estimativa de condicionamento

Chama `mat4_cond_estimate` (razão max/min de pivôs por eliminação gaussiana). Quando singular, retorna 0.0 (sentinel inicializado em `out->cond_estimate = 0.0`). Correto.

### 1.4. Correção aplicada nesta auditoria

**Bug encontrado e corrigido em `apps/reproduce_cube_sparse.c`:**

O código original avaliava o diagnóstico em pontos incorretos:

```c
// ERRADO (antes da correção)
double hx = 0.5 * L / (double)nx_cells;
mls_connectivity_stats_uniform_grid(nodes, node_count,
                                    hx, L - hx, nx_cells, ...);
```

`mls_connectivity_stats_uniform_grid` com `[xmin, xmax, n]` coloca pontos em
`xmin + (ix + 0.5) * (xmax - xmin) / n`. Usando `xmin = hx = L/(2*nx_cells)`:

```
px_ix0 = 0.644  (deveria ser 0.333 para nx_cells=15, L=10)
```

O desvio no primeiro ponto é ~L/(2*nx_cells) ≈ 0.33 para nx_cells=15, que é
metade de uma célula. Os pontos estavam sistematicamente deslocados para o interior
do domínio, inflando artificialmente o número de nós ativos e escondendo casos de
baixa conectividade próximos à borda.

**Correção aplicada:**

```c
// CORRETO (após a correção)
mls_connectivity_stats_uniform_grid(nodes, node_count,
                                    0.0, L, nx_cells,
                                    0.0, L, ny_cells,
                                    0.0, L, nz_cells,
                                    &diag);
```

Resultado antes da correção para 11x11x11:

```
min active nodes: 27  (artificialmente alto)
```

Resultado após a correção para 11x11x11:

```
min active nodes:  8  (valor real — células de borda têm menos nós ativos)
mean active nodes: 23.56
max cond estimate: 6.76e+02
```

Conclusão: todas as células têm active_nodes ≥ 8 (acima do limiar de 4). Nenhuma
célula inválida. O pipeline pode avançar.

---

## 2. Base MLS linear 3D — VERIFICADO CORRETO

Base implementada: `p^T(x) = [1, x, y, z]` — conforme artigo, Eq. (8)-(9).

O solver local `mat4_solve` usa eliminação gaussiana com pivotamento parcial.
Tolerância de pivô: `1e-12`.

Gradientes implementados em `src/mls.c`:

```
dPhi_I/dx = dp^T/dx * A^{-1}B_I + p^T * d(A^{-1}B_I)/dx
           = q_I[1] + p^T(x) * r_dx_I
```

onde `q_I = A^{-1}B_I` e `r_dx_I = A^{-1}(dB_I/dx - dA/dx * q_I)`. Implementação
verificada linha a linha contra a fórmula. Correto.

---

## 3. Sistema aumentado [K G^T; G 0] — VERIFICADO CORRETO

### 3.1. Versão densa (`lagrange_assemble_augmented_system_dense`)

```c
/* K block: A_aug[row, col] = K[row, col]  for row,col in 0..n-1  */
A_aug[(row * total_size) + col] = K[(row * node_count) + col];

/* G^T block: A_aug[node, node_count+constraint] = G[constraint, node]  */
A_aug[(node * total_size) + aug_row] = value;

/* G block: A_aug[node_count+constraint, node] = G[constraint, node]  */
A_aug[(aug_row * total_size) + node] = value;
```

Posicionamento correto.

### 3.2. Versão esparsa (`lagrange_assemble_augmented_sparse_coo`)

```c
/* G^T: A_aug[n, node_count+c] = G[c, n]  ✓ */
sparse_coo_add(A_aug_coo, n, node_count + c, val);

/* G:   A_aug[node_count+c, n] = G[c, n]  ✓ */
sparse_coo_add(A_aug_coo, node_count + c, n, val);
```

Posicionamento correto. **RISCO DE DESEMPENHO: esta função insere TODAS as entradas
de G, incluindo zeros** (detalhado na seção de riscos abaixo).

### 3.3. Bloco G: definição de G[c, n]

`G[c, n] = Phi_n(x_c)` onde `x_c` é o ponto Dirichlet c.

Montado em `src/dirichlet.c`:

```c
G[(point_index * node_count) + node_index] = shape_values[value_index].phi;
```

Correto. Apenas nós com suporte ativo em `x_c` recebem entradas não-nulas.

### 3.4. Vetor do lado direito

```
b_aug = [F; q]  onde  F_i = 0 (sem fonte volumétrica, rho=0)
                      q_c = V_d(x_c)  (potencial prescrito)
```

Implementado corretamente. `F` é zero (calloc). `q` recebe `point->value` de
cada `DirichletPoint`.

### 3.5. Condições de contorno do cubo

Implementado em `cube_generate_dirichlet_points_from_nodes`:

- `z = L` → `value = V0 = 10`  (face superior)
- todas as outras faces → `value = 0`
- face superior tem prioridade sobre arestas superiores

Correto segundo `docs/06_condicoes_de_contorno_e_lagrange.md`.

---

## 4. GMRES — VERIFICADO CORRETO

### 4.1. Interface

```c
int gmres_solve(const SparseCSR *A, const double *b, double *x,
                double tol, int max_iter, int restart, GmresResult *result);
```

Recebe CSR. Retorna `GmresResult`:

```c
typedef struct GmresResult {
    int    converged;       /* 1/0 */
    int    iterations;      /* total Arnoldi steps (todos os restarts) */
    double residual_init;   /* ||b - A*x0||  */
    double residual_final;  /* ||b - A*x_k||  (recomputado ao final) */
} GmresResult;
```

Todas as informações necessárias estão presentes.

### 4.2. Critério de convergência

Tolerância relativa: convergência quando `||r_k|| / ||r_0|| < tol`.
Verificação dentro do ciclo de Arnoldi: `fabs(g[j+1]) / beta0 < tol`.
Verificação no restart: `beta / beta0 < tol`.

Resíduo final reportado é calculado explicitamente como `||b - A*x||` após a saída.
Correto.

### 4.3. Resultados medidos para 11x11x11

```
restart=200, max_iter=10000, tol=1e-9
iterations=67, converged=YES
residual_init=1.10e+02, residual_final=8.70e-08
rel_residual=7.91e-10
solution_time=0.008s
```

### 4.4. Risco GMRES para casos maiores

Para 15x15x15: augmented_size ≈ 4553. O restart atual de 200 pode ser insuficiente
e `max_iter=10000` pode não ser alcançado. Ver recomendações.

---

## 5. Solução analítica e comparação de erro — CORRETO COM RESSALVAS

### 5.1. Fórmula analítica

`src/analytical.c` implementa Eqs. (13)-(15) do artigo:

```c
V = Σ_n Σ_m a_mn * sin(n*pi*x/L) * sin(m*pi*y/L) * sinh(k_mn*z)
```

com `n, m` ímpares. Correto.

### 5.2. Métrica de erro implementada

A aplicação usa norma L2 discreta em grade de amostragem:

```
e_global  = sqrt( Σ (v_num - v_exact)^2 / Σ v_exact^2 )
e_interior = idem, excluindo pontos de contorno
```

O artigo usa a integral contínua (Eq. 16):

```
e = sqrt( ∫(V_EFG - V_exact)^2 dΩ / ∫ V_exact^2 dΩ )
```

A métrica discreta é uma aproximação razoável. Não são comparáveis numericamente
com os resultados do artigo, mas a tendência de convergência é válida.

### 5.3. Max abs error ≈ V0 é ESPERADO e não é bug

O erro máximo absoluto ≈ 10 = V0 ocorre nos pontos de aresta superior do cubo,
como (0, y, L) ou (L, y, L), onde as condições de Dirichlet são contraditórias:

- Face `z=L` impõe `V = V0`
- Faces `x=0` e `x=L` impõem `V = 0`

A série analítica resolve essa contradição via continuidade: em (0, y, L),
`sin(n*pi*0/L) = 0` para todo n, portanto a série retorna V_exact = 0.
O EFG com multiplicadores de Lagrange impõe `V = V0` na face superior COM
PRIORIDADE, de modo que o ponto (0, y, L) recebe `V_d = V0` na montagem.

Esta discrepância é intrínseca ao benchmark e não representa erro de implementação.
**O indicador relevante é o erro relativo interior**, que mede apenas pontos
estritamente dentro do cubo.

### 5.4. Resultados medidos

| Caso | Interior rel. error | Artigo (não-uniforme) |
|------|--------------------|-----------------------|
| 5x5x5 uniform  | 21.8% | — |
| 11x11x11 uniform | 5.5% | 1.40% |

A diferença entre 5.5% e 1.40% é esperada:
- O artigo usa nuvem não-uniforme com mais nós nas arestas superiores
- O artigo mede sobre o grid não-uniforme, não em amostragem discreta regular

---

## 6. Riscos técnicos

### Risco 1 — MODERADO: `lagrange_assemble_augmented_sparse_coo` insere zeros

A função itera sobre TODOS os `node_count * constraint_count` pares (n, c) e
insere as entradas em COO mesmo quando `G[c, n] = 0`. Para 15x15x15:

```
2 * 3375 * 1178 = 7.950.900 entradas inseridas no COO
Mas apenas ~3375 * 1178 / 1331 * 16 ≈ ~48.000 são não-zero
Overhead: ~165x
```

O CSR resultante contém zeros explícitos que tornam o matvec do GMRES mais lento.

**Impacto:** Esta função NÃO é usada por `reproduce_cube_sparse.c` para casos
grandes (o app tem sua própria montagem com filtro de zeros). Ela é usada apenas
nos testes `test_lagrange_augmented_sparse` e `test_gmres` (caso 3x3x3). Para esses
casos pequenos, o overhead é aceitável.

**Ação:** Não corrigir agora. Documentar e evitar usar para casos maiores.

### Risco 2 — MODERADO: estimativa de capacidade inicial do K_coo subestimada

```c
int k_cap = node_count * 30 + 16;  // em reproduce_cube_sparse.c
```

Para 11x11x11: k_cap = 39.946, K nnz real = 117.649 (3× maior).
Para 15x15x15: k_cap ≈ 101.266, K nnz real estimado ≈ 297.000 (3× maior).

`sparse_coo_add` dobra capacidade automaticamente, então não é bug. É ineficiência:
3+ realocações em vez de 1 por caso grande.

**Ação recomendada ao Codex:** Ajustar estimativa:

```c
// Para grade uniforme NxNxN, K é praticamente denso dentro do suporte:
// cada nó tem ~27 vizinhos no suporte → ~27^2 ≈ 729 entradas por nó,
// mas com compressão de duplicatas ficam ~88 entradas/nó (medido em 11x11x11).
// Usar 100 * node_count como estimativa conservadora.
int k_cap = node_count * 100 + 16;
```

### Risco 3 — ALTO PARA CASOS MUITO GRANDES: teto de memória por K denso

O pipeline atual monta K como matriz densa antes de converter para COO:

```
K_dense = calloc(node_count^2, sizeof(double))
```

| Grade | Nós | K denso | G denso | Total |
|-------|-----|---------|---------|-------|
| 11x11x11 | 1.331 | 14 MB | 6.4 MB | 21 MB |
| 13x13x13 | 2.197 | 38 MB | 14 MB | 52 MB |
| 15x15x15 | 3.375 | 91 MB | 31 MB | 122 MB |
| 17x17x17 | 4.913 | 193 MB | 58 MB | 251 MB |
| 21x21x21 | 9.261 | 686 MB | 146 MB | 832 MB |

Para 15x15x15, ~122 MB é viável em sistemas com >512 MB RAM. Acima de 17x17x17,
o pipeline denso se torna problemático.

**Ação:** Aceitar para agora (15x15x15 é o alvo). Para casos maiores no futuro,
o Codex deve implementar montagem direta K esparso sem passar por K denso.

### Risco 4 — BAIXO: condicionamento do sistema cresce com o refinamento

```
11x11x11: max_cond_estimate = 6.76e+02
```

O condicionamento da matriz de momento A(x) cresce com a razão de escalas. Para
malhas mais finas com o mesmo d_mI = 1.5 × 3o_vizinho, as dimensões absolutas
dos suportes diminuem mas o ratio de escalas (d_mI na base p^T = [1, x, y, z])
cresce. Isso pode causar mal condicionamento da matriz de momento para grids muito finos.

GMRES com tolerância 1e-9 e sistema aumentado simétrico-indefinido pode degradar
para casos muito mal condicionados. Monitorar o campo `max_cond_estimate` nos relatórios.

**Critério de bloqueio:** Se `max_cond_estimate > 1e6` em algum caso, investigar
antes de prosseguir.

### Risco 5 — BAIXO: parâmetros GMRES para casos maiores

O caso 11x11x11 converge em 67 iterações com `restart=200`. Para 15x15x15
(sistema augmentado ≈ 4553), é possível que sejam necessárias mais iterações.

**Parâmetros recomendados para 15x15x15:**

```c
gmres_tol     = 1e-9
gmres_restart = 300
gmres_max_iter = 20000
```

Se não convergir: aumentar `restart` antes de aumentar `max_iter` (Krylov maior
ajuda mais por iteração do que mais restarts com Krylov pequeno).

---

## 7. Critérios de aceitação para o próximo caso

Um caso é considerado PASSOU se satisfaz TODOS os seguintes critérios:

```
support_lt_4        = 0          (nenhuma célula inválida)
support_lt_8        = 0          (preferível; se > 0, investigar localização)
mls_failures        = 0          (nenhuma falha de matriz de momento)
max_cond_estimate  < 1e6         (sistema momento não degenerado)
gmres_converged     = YES
rel_residual       < 1e-8        (tolerância efetiva)
rel_error_interior reportado
```

O critério de erro interior não é bloqueio nesta fase (o objetivo é verificar
o pipeline, não reproduzir o artigo ainda).

---

## 8. Checklist para o Codex

### 8.1. Ajustes obrigatórios antes de rodar casos maiores

- [x] Corrigir posicionamento dos pontos de diagnóstico em `reproduce_cube_sparse.c`
      (feito nesta auditoria — commit pendente)

- [ ] Ajustar estimativa de capacidade K_coo:
      Trocar `node_count * 30` por `node_count * 100`

- [ ] Ajustar parâmetros GMRES para casos maiores:
      `restart=300, max_iter=20000` para ≥ 13x13x13

- [ ] Adicionar caso 13x13x13 e 15x15x15 em `main()` de `reproduce_cube_sparse.c`
      com `run_dense_cmp=0` e parâmetros GMRES ajustados

### 8.2. Comandos que o Codex deve executar

```bash
# Passo 1: build completo
cmake --build build

# Passo 2: suite de testes
ctest --test-dir build --output-on-failure

# Passo 3: rodar todos os casos
./build/reproduce_cube_sparse

# Passo 4: registrar saída em docs/dev/TEST_REPORT.md
```

### 8.3. O que registrar em TEST_REPORT.md para cada caso

```text
nodes:
constraints:
augmented_size:
K_nnz:
A_aug_nnz:
support_min:
support_mean:
support_max:
support_lt_4:       # DEVE ser 0 para prosseguir
support_lt_8:
max_cond_estimate:
mls_failures:       # DEVE ser 0 para prosseguir
gmres_converged:    # DEVE ser YES para prosseguir
gmres_iterations:
residual_init:
residual_final:
rel_residual:
rel_error_global:
rel_error_interior:
max_abs_error:
assembly_time_s:
solve_time_s:
```

---

## 9. O que NÃO deve ser feito ainda

1. **Nuvens não uniformes** — pausado conforme TODO.md seção 11.
2. **Precondicionador GMRES** — não implementar antes de confirmar convergência sem ele.
3. **Montagem G esparsa nativa** — não é necessário para 15x15x15; traria complexidade.
4. **Octree ou k-d tree** — O(n^2) de `support_assign_article_default` é aceitável para n ≤ 3375.
5. **Refatoração de APIs** — todas as interfaces estão corretas; não mexer.
6. **Figura do artigo (plano x=5.33)** — Implementar APENAS após confirmar resultados nos casos maiores.
7. **Exportação CSV de diagnóstico por ponto** — mencionada no TODO mas desnecessária agora.

---

## 10. Sequência recomendada de execução

```
1. Corrigir k_cap → k_cap = node_count * 100 + 16
2. Adicionar casos 13x13x13 e 15x15x15 em reproduce_cube_sparse.c
3. Ajustar GMRES restart=300, max_iter=20000 para esses casos
4. cmake --build build && ctest
5. ./build/reproduce_cube_sparse  →  registrar em TEST_REPORT.md
6. Verificar critérios de aceitação
7. Se 15x15x15 passar → implementar exportação do plano x=5.33
```

---

## 11. Recomendação final

```
AVANÇAR COM AJUSTES PONTUAIS
```

A base está correta. O pipeline completo (sparse K + sparse aug + GMRES) funciona
para 11x11x11 e converge em 67 iterações com resíduo relativo 7.9e-10. Os ajustes
necessários são pequenos (k_cap, parâmetros GMRES, novos casos) e não envolvem
risco de regressão nos testes existentes.

O único bug real encontrado (pontos de diagnóstico deslocados) foi corrigido nesta
auditoria. Todos os 28 testes continuam passando após a correção.
