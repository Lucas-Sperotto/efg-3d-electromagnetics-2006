# TODO — Fechamento científico e computacional do EFG 3D

Este arquivo organiza o estado atual e o plano de continuidade do projeto `efg-3d-electromagnetics-2006`.

O objetivo principal é reproduzir, de forma didática e verificável, o problema eletrostático tridimensional do artigo:

> The Element-Free Galerkin Method in Three-Dimensional Electromagnetic Problems  
> IEEE Transactions on Magnetics, 2006.

A prioridade agora não é criar investigações secundárias, mas consolidar a formulação numérica para permitir refinamento suficiente e reprodução confiável das figuras do artigo.

---

## 0. Estado atual do projeto

### 0.1. Base já implementada

- [x] Estrutura C11 com CMake.
- [x] Funções peso escalares.
- [x] Derivadas escalares das funções peso.
- [x] Peso tensorial 3D com gradiente físico.
- [x] Cálculo de suporte nodal baseado no terceiro vizinho multiplicado por 1,5.
- [x] Solver local 4x4 para MLS linear 3D.
- [x] MLS linear 3D sem derivadas.
- [x] MLS linear 3D com gradientes.
- [x] Quadratura de Gauss 2x2x2 como padrão.
- [x] Regras Gauss-Legendre 1D/3D para ordens 1..8 em estudos controlados.
- [x] Montagem local da matriz de rigidez.
- [x] Montagem global densa da matriz de rigidez.
- [x] Montagem densa das restrições de Dirichlet via multiplicadores de Lagrange.
- [x] Montagem densa do sistema aumentado:

```text
  [ K  G^T ]
  [ G   0  ]
```

- [x] Solver denso didático por eliminação de Gauss.
- [x] Teste integrado artificial com solução linear.
- [x] Gerador de problema cúbico regular.
- [x] Geração de Dirichlet a partir de grade regular.
- [x] Geração de Dirichlet a partir de nuvem arbitrária de nós.
- [x] Aplicação `reproduce_cube_smoke`.
- [x] Script de análise de erro por região.
- [x] Script de estudo de refinamento regular.
- [x] Script de plotagem da convergência.
- [x] Primeira tentativa de nuvem não uniforme refinada no topo.
- [x] Diagnóstico de conectividade MLS (`mls_diagnostic`).
- [x] Estimativa de condicionamento da matriz de momento A(x).
- [x] Módulo de matriz esparsa COO/CSR (`sparse_matrix`).
- [x] Montagem esparsa da rigidez (`global_stiffness_assemble_sparse_coo`).
- [x] Montagem esparsa do sistema aumentado (`lagrange_assemble_augmented_sparse_coo`).
- [x] Solver iterativo GMRES reiniciado (`gmres_solve`).
- [x] Solver LAPACK denso opcional (`solve_augmented_lapack_dense`) para validar
  o sistema aumentado esparso contra `dgesv`.
- [x] Norma relativa integral da Eq. (16) no domínio (`relative_error_norm_domain_integral`).
- [x] Estudo de sensibilidade da ordem de integração para `refine15`.

### 0.2. Evidências recentes (29/29 testes passando)

- `test_mls_diagnostic` — conectividade e condicionamento em grade 3×3×3 e 5×5×5.
- `test_sparse_matrix` — inserção, duplicatas, compressão, matvec CSR.
- `test_sparse_global_stiffness` — sparse K ≈ dense K para 3×3×3 (tol 1e-10).
- `test_lagrange_augmented_sparse` — sparse A_aug ≈ dense A_aug para 5×5 manual e 3×3×3 cube.
- `test_gmres` — solução 5×5 conhecida + GMRES vs. solver denso no cubo 3×3×3.
- `test_lapack_dense_solver` — conversão CSR -> dense column-major e solve LAPACK 5×5 quando disponível.
- `test_error_norm` — norma discreta e norma integral Eq. (16) com campos constantes,
  reprodução linear, denominador zero e falha MLS.
- `test_gauss` — regras Gauss-Legendre 1..8, soma dos pesos 1D/3D,
  integração de volume físico, função linear e polinômio quadrático.

---

## 1. Diagnóstico importante

A primeira tentativa de nuvem não uniforme não melhorou o resultado.

Comparação observada:

```text
Regular 11x11x11 com integração 15x15x15:
nodes: 1331
constraints: 602
internal relative error: ~0.05294

Não uniforme base_n=7 refine_n=11:
nodes: 1308
constraints: 786
internal relative error: ~0.15300
```

Conclusão provisória:

- a nuvem não uniforme ficou dominada por pontos de contorno;
- houve aumento relativo das restrições de Dirichlet;
- a distribuição pode ter piorado o condicionamento local;
- ainda não há diagnóstico de cardinalidade e condicionamento do MLS;
- ainda não há estrutura esparsa para resolver casos suficientemente refinados.

Portanto, a linha de trabalho com nuvens não uniformes deve ficar pausada.

---

## 2. Correção de rumo

Antes de continuar refinando a geometria ou tentando reproduzir a figura final, é necessário garantir:

1. cada ponto de integração possui número suficiente de nós ativos no suporte;
2. a cardinalidade local não é inferior ao número de termos da base;
3. a matriz de momento MLS não fica singular ou mal condicionada;
4. a montagem global deve migrar de densa para esparsa;
5. o sistema aumentado deve ser resolvido por método iterativo adequado, inicialmente GMRES;
6. só depois disso deve-se buscar refinamento suficiente para reproduzir a figura do artigo.

Para a base linear 3D:

```text
p(x,y,z) = [1, x, y, z]
```

a cardinalidade mínima teórica é:

```text
4 nós ativos
```

mas, na prática, é recomendável ter mais do que 4 nós no suporte para reduzir singularidade e mau condicionamento.

---

## 3. Prioridade 1 — Diagnóstico de conectividade MLS

### 3.1. Objetivo

Criar ferramentas para avaliar a conectividade local de cada ponto de integração e de cada ponto de avaliação.

Para cada ponto analisado, calcular:

- número de nós ativos no suporte;
- se a cardinalidade é menor que 4;
- se a cardinalidade é pequena, por exemplo menor que 8;
- se a matriz de momento MLS é singular ou quase singular;
- distribuição estatística da conectividade.

### 3.2. Tarefas

- [x] Criar módulo de diagnóstico MLS.
- [x] Criar função para contar nós ativos no suporte de um ponto.
- [x] Criar função para avaliar estatísticas de conectividade em uma grade de pontos.
- [x] Criar teste unitário para conectividade em grade regular.
- [ ] Criar teste unitário para conectividade em nuvem não uniforme.
- [ ] Criar app ou script que gere CSV com:

  ```text
  x,y,z,active_nodes,status
  ```

- [ ] Gerar relatório para os casos:

  ```text
  3x3x3
  5x5x5
  7x7x7
  9x9x9
  11x11x11
  ```

### 3.3. Critério de parada

Não avançar para refinamentos maiores enquanto houver pontos de integração ou amostragem com:

```text
active_nodes < 4
```

Também registrar pontos com:

```text
active_nodes < 8
```

como região numericamente suspeita.

---

## 4. Prioridade 2 — Diagnóstico da matriz de momento MLS

### 4.1. Objetivo

Avaliar a qualidade da matriz:

```text
A(x) = Σ_I w_I p_I p_I^T
```

em pontos de integração e amostragem.

### 4.2. Tarefas

- [x] Criar função para montar `A(x)` sem calcular necessariamente as funções de forma.
- [x] Criar função para estimar condicionamento ou pelo menos detectar pivôs pequenos durante a eliminação.
- [ ] Registrar:

  - menor pivô;
  - falha de solução;
  - número de nós ativos;
  - coordenada do ponto;
  - tipo de ponto: integração, avaliação ou Dirichlet.
- [ ] Gerar CSV diagnóstico.
- [ ] Comparar grade regular e nuvem não uniforme.

### 4.3. Critério de parada

Não buscar reprodução quantitativa enquanto houver falhas frequentes em MLS ou regiões com matriz de momento quase singular.

---

## 5. Prioridade 3 — Matriz esparsa

### 5.1. Justificativa

A montagem densa permitiu validar o pipeline, mas não permite refinamento suficiente.

Exemplo:

```text
15x15x15 nós = 3375 nós
restrições de superfície ≈ 1178
sistema aumentado ≈ 4553 x 4553
```

Uma matriz densa desse porte é cara em memória e muito cara para eliminação de Gauss.

O artigo usa GMRES porque o sistema aumentado com multiplicadores de Lagrange é grande e não é positivo definido.

### 5.2. Estrutura esparsa mínima

Implementar primeiro formato COO didático:

```c
typedef struct {
    int row;
    int col;
    double value;
} SparseEntry;
```

Depois converter para CSR:

```text
row_ptr
col_ind
values
```

### 5.3. Tarefas

- [x] Criar módulo `sparse_matrix`.
- [x] Implementar estrutura COO.
- [x] Permitir acumular entradas repetidas.
- [x] Implementar ordenação por `(row, col)`.
- [x] Implementar compressão de duplicatas.
- [x] Implementar conversão COO → CSR.
- [x] Implementar produto matriz-vetor CSR.
- [x] Criar testes unitários para:

  - inserção;
  - duplicatas;
  - compressão;
  - produto matriz-vetor;
  - comparação com matriz densa pequena.

---

## 6. Prioridade 4 — Montagem esparsa da rigidez

### 6.1. Objetivo

Criar equivalente esparso de:

```c
global_stiffness_assemble_dense(...)
```

sem remover a versão densa.

### 6.2. Tarefas

- [x] Criar `global_stiffness_assemble_sparse_coo(...)`.
- [x] Reutilizar a montagem local já existente.
- [x] Acumular entradas locais em COO.
- [x] Comprimir duplicatas.
- [x] Converter para CSR.
- [x] Comparar, em problemas pequenos, a matriz esparsa com a matriz densa.
- [x] Criar teste:

  ```text
  sparse K ≈ dense K
  ```

para `3x3x3`.

---

## 7. Prioridade 5 — Sistema aumentado esparso

### 7.1. Objetivo

Montar em formato esparso:

```text
[ K  G^T ]
[ G   0  ]
```

### 7.2. Tarefas

- [x] Criar montagem esparsa de `G`.
- [x] Criar montagem esparsa do sistema aumentado.
- [x] Comparar com o sistema denso em caso pequeno.
- [x] Criar teste:

  ```text
  sparse A_aug * x ≈ dense A_aug * x
  ```

---

## 8. Prioridade 6 — GMRES

### 8.1. Objetivo

Implementar GMRES para resolver:

```text
A_aug x = b_aug
```

com matriz esparsa CSR.

### 8.2. Tarefas

- [x] Implementar GMRES básico sem precondicionador.
- [x] Criar teste com sistema pequeno conhecido.
- [x] Testar contra `dense_solve(...)` em sistema pequeno.
- [x] Registrar:

  - número de iterações;
  - resíduo inicial;
  - resíduo final;
  - tolerância;
  - convergência ou falha.
- [x] Implementar versão reiniciada, por exemplo GMRES(m), se necessário.
- [x] Avaliar necessidade de precondicionador.

### 8.3. Possíveis precondicionadores futuros

- [x] Jacobi — implementado; **não eficaz** para o sistema aumentado [K G^T; G 0].
  O bloco (2,2) tem diagonal zero → diag_inv = 1 para as linhas dos multiplicadores →
  disparidade de escala não melhora o condicionamento. 96 iter (vs. 69 sem precond)
  para regular 15×15×15; 858 iter (vs. 832) para nonuniform_refine. API mantida.
- [ ] Precondicionador bloco-diagonal [diag(K); G diag(K)⁻¹ G^T] — abordagem correta
  para sistemas de ponto-sela; requer montagem explícita do complemento de Schur aproximado.
- [ ] SSOR ou escalonamento simétrico — pode ajudar sem bloco-diagonal.
- [ ] ILU simples — mais caro, mas mais geral.

### 8.4. Validação direta opcional com LAPACK

- [x] Criar opção de solver:

  ```bash
  ./build/reproduce_cube_sparse --case refine15 --solver gmres
  ./build/reproduce_cube_sparse --case refine15 --solver lapack-dense
  ```

- [x] Implementar conversão CSR -> matriz densa column-major.
- [x] Resolver com LAPACK `dgesv` quando disponível.
- [x] Manter build funcional sem LAPACK.
- [x] Registrar tempo de conversão, tempo de fatoração/solução, resíduo final
  e diferença relativa GMRES/LAPACK.
- [x] Validar `refine15`: diferença relativa GMRES/LAPACK `8.285314e-10`.
- [x] Validar `nonuniform_refine`: diferença relativa GMRES/LAPACK `3.368672e-08`.

### 8.5. Norma relativa integral da Eq. (16)

- [x] Implementar integral de domínio:

  ```text
  e = sqrt( int_Omega (V_EFG - V_exato)^2 dOmega
            / int_Omega V_exato^2 dOmega )
  ```

- [x] Avaliar `V_EFG = Σ_I Phi_I(x_gp) u_I` nos pontos de Gauss.
- [x] Usar Gauss 2x2x2 como padrão para manter compatibilidade com a montagem.
- [x] Permitir ordem de norma explícita para separar avaliação de erro da montagem.
- [x] Tornar `relative_error_global` a métrica principal da Eq. (16).
- [x] Preservar a métrica antiga como `relative_error_discrete_global`.
- [x] Validar `refine15`: Eq. (16) `8.217121e-02`.
- [x] Validar `nonuniform_refine`: Eq. (16) `8.035461e-02`.

### 8.6. Estudo controlado da ordem de integração

- [x] Generalizar a quadratura Gauss-Legendre para ordens 1..8.
- [x] Manter ordem 2 como padrão de montagem e norma histórica.
- [x] Separar ordem da montagem de K e ordem da norma Eq. (16).
- [x] Gerar `data/output/integration_order_norm_sensitivity_refine15.csv`.
- [x] Gerar `data/output/integration_order_solution_sensitivity_refine15.csv`.
- [x] Validar estudo de norma: ordens 5..8 estabilizam em torno de `9.18e-02`.
- [x] Validar estudo de montagem: ordens 2..5 estabilizam em torno de `9.20e-02`
  usando norma de referência ordem 6.
- [x] Registrar que montagem ordem 1 converge, mas piora o erro integral para
  `9.842689e-02`.

### 8.7. Análise regional da norma integral da Eq. (16)

- [x] Criar CLI:

  ```bash
  ./build/reproduce_cube_sparse --case refine15 --error-region-study
  ```

- [x] Fixar `refine15`: nós 15x15x15, células 15x15x15, montagem ordem 2,
  norma ordem 6 e GMRES.
- [x] Gerar `data/output/error_region_integral_refine15.csv`.
- [x] Avaliar regiões: domínio completo, interior aberto, núcleos com delta,
  caixa central, camada superior, bandas de aresta superiores e interior da
  face superior.
- [x] Validar `mls_failures = 0` em todas as regiões.
- [x] Registrar região dominante: `upper_layer` concentra 98,24 % do numerador
  global; `upper_edge_bands` concentra 97,14 % do numerador global.
- [x] Registrar que o núcleo do domínio não domina o erro: `central_box`
  contribui apenas 0,14 % do numerador global.

---

## 9. Prioridade 7 — Reprodução do cubo com refinamento suficiente

### 9.1. Objetivo

Resolver o problema do cubo com refinamento suficiente para aproximar a figura do artigo.

### 9.2. Casos-alvo

Casos uniformes iniciais:

- [x] 11x11x11 nós, 15x15x15 células — executado com `reproduce_cube_sparse`.
- [x] 13x13x13 nós, 15x15x15 células — executado com `reproduce_cube_sparse --case refine13`.
- [x] 15x15x15 nós, 15x15x15 células — executado com `reproduce_cube_sparse --case refine15`.

Nuvem não uniforme:

- [x] base_n=11, top_n=13, n_extra_slices=4, z_frac=0.30 — executado com `reproduce_cube_sparse --case nonuniform`.
  - 1974 nós, 762 restrições (38,6 % de superfície), erro interior 2,96 %, GMRES 391 iterações.

### 9.3. Métricas

- [x] erro relativo discreto global;
- [x] erro relativo integral de domínio pela Eq. (16);
- [x] erro relativo discreto interno;
- [x] erro máximo absoluto;
- [x] erro no plano `x = 5.33`;
- [x] erro em regiões próximas às arestas superiores;
- [x] número de iterações GMRES;
- [x] tempo de montagem;
- [x] tempo de solução;
- [x] tempo de conversão CSR -> dense e tempo `dgesv` para validação LAPACK;
- [x] diferença relativa GMRES/LAPACK nos casos regular e não uniforme refinado;
- [ ] memória aproximada.

---

## 10. Prioridade 8 — Reprodução da figura

### 10.1. Figura de interesse

A figura principal deve ser reproduzida em plano semelhante ao usado no artigo, especialmente:

```text
x = 5.33
```

### 10.2. Tarefas

- [x] Criar modo no app para exportar plano `x = 5.33`.
- [x] Gerar CSV:

  ```text
  y,z,V_num,V_exact,abs_error
  ```

  Concluído em `data/output/cube_plane_x_5_33_refine15.csv` com colunas:

  ```text
  x,y,z,V_num,V_exact,abs_error
  ```

- [x] Criar script Python para plotar:

  - solução numérica;
  - solução analítica;
  - erro absoluto;
  - mapa de contorno.

  Concluído em `scripts/plot_cube_plane.py`, gerando figuras e métricas em:

  ```text
  data/output/figures/
  ```

- [x] Ajustar script de contorno para comparação com a Figura 3:

  - iso-linhas de `V_num`;
  - iso-linhas de `V_exact`;
  - sobreposição com `V_num` sólido e `V_exact` tracejado;
  - figura estilo artigo com fundo branco e níveis fixos `1,2,3,4,5,6,7,8,9,10`.

- [x] Comparar erro máximo com valor reportado no artigo.
- [x] Documentar divergências.

  Concluído em `docs/figure3_reproduction.md`: a figura principal fica com o
  caso regular `15x15x15`; a nuvem `nonuniform` é comparação suplementar.
  Os valores do artigo (`1,40 %` e `0,15 V`) ainda não foram reproduzidos
  quantitativamente. A divergência foi atribuída à grade/nuvem ainda diferente
  da publicação, série truncada e ausência de comparação digital pixel a pixel.

---

## 11. Prioridade 9 — Nuvens não uniformes

### 11.1. Status

A primeira tentativa de nuvem não uniforme (base_n=7, refine_n=11) foi diagnosticada como
falha por excesso de nós de superfície (60 % de restrições). Uma segunda abordagem foi
implementada e validada com `cube_generate_article_cloud`.

### 11.2. Resultados da nuvem article cloud (2026-05-06)

| Parâmetro | Valor |
| --- | --- |
| base_n | 11 |
| top_n | 13 |
| n_extra_slices | 4 |
| z_frac | 0,30 |
| Nós totais | 1974 |
| Restrições | 762 (38,6 % de superfície) |
| support_lt_4 | 0 |
| mls_failures | 0 |
| GMRES convergiu | YES |
| GMRES iterações | 391 |
| Erro relativo interior | 2,96 % |
| Erro relativo interior 11×11×11 regular | 5,13 % |
| Erro relativo interior 15×15×15 regular | 2,80 % |

A nuvem reduziu o erro interior em relação ao 11×11×11 regular com número semelhante de nós
e ficou próxima do 15×15×15 regular com menos nós. A estratégia de adicionar somente nós
interiores-em-xy na zona superior evita o problema de dominância por restrições.

### 11.3. Refinamento da nuvem article cloud (2026-05-06)

| Parâmetro | Valor |
| --- | --- |
| base_n | 13 |
| top_n | 15 |
| n_extra_slices | 4 |
| z_frac | 0,30 |
| Nós totais | 3089 |
| Restrições | 1082 |
| support_lt_4 | 0 |
| mls_failures | 0 |
| GMRES convergiu | YES |
| GMRES iterações | 832 |
| Erro relativo interior 3D | 2,78 % |
| Erro relativo interior do plano | 4,44 % |
| Erro relativo na janela central | 2,90 % |

O refinamento melhorou levemente o erro 3D interior em relação ao 15×15×15 regular
(`2,78 %` contra `2,80 %`) e reduziu o erro relativo interior do plano (`4,44 %`
contra `4,55 %`). Porém, piorou a janela central e elevou muito o custo do GMRES
(`832` iterações). Portanto, ainda não substitui o caso regular como figura
principal.

### 11.4. Próximos passos opcionais

- [ ] Ajustar z_frac e n_extra_slices para concentrar mais nós na região de gradiente crítico.
- [x] Testar base_n=13, top_n=15 para avaliar convergência da nuvem não uniforme.
- [ ] Avaliar necessidade de precondicionador para reduzir as 832 iterações do GMRES.
- [ ] Evitar nuvens dominadas por pontos de superfície (razão documentada no §1).

---

## 12. O que estava sendo proposto antes, mas deve ficar pausado

As seguintes ideias estavam no caminho, mas devem ser pausadas:

- [ ] Comparar estatisticamente CSV regular versus não uniforme.
- [ ] Ajustar empiricamente a nuvem não uniforme.
- [ ] Criar muitos gráficos de erro por região.
- [ ] Refinar heurísticas de concentração de nós no topo.
- [ ] Investigar overshoot em nuvens não uniformes.

Essas tarefas são úteis, mas são secundárias. O gargalo principal agora é estrutural:

```text
conectividade MLS
condicionamento
matriz esparsa
GMRES
refinamento suficiente
```

---

## 13. Ordem recomendada de execução a partir de agora

1. Diagnóstico de conectividade MLS.
2. Diagnóstico da matriz de momento MLS.
3. Matriz esparsa COO/CSR.
4. Montagem esparsa de K.
5. Montagem esparsa de G.
6. Sistema aumentado esparso.
7. GMRES básico.
8. Comparação GMRES versus solver denso em casos pequenos.
9. Rodar casos uniformes maiores.
10. Exportar plano `x = 5.33`.
11. Reproduzir figura do artigo.
12. Só depois voltar para nuvens não uniformes.

---

## 14. Critério de sucesso intermediário

Antes de buscar a figura final, o projeto deve conseguir rodar:

```text
15x15x15 nós
15x15x15 células de integração
Dirichlet via Lagrange
GMRES
exportação de plano x = 5.33
comparação com solução analítica
```

com relatório contendo:

```text
nodes
constraints
nonzeros
GMRES iterations
residual
relative error
max abs error
runtime
```

---

## 15. Critério de sucesso final

O projeto será considerado cientificamente consistente quando conseguir:

- reproduzir qualitativamente a distribuição de potencial no cubo;
- obter erro interno decrescente com refinamento;
- exportar figura comparável à do artigo;
- documentar diferenças em relação ao artigo;
- explicar claramente o papel de:

  - suporte nodal;
  - cardinalidade local;
  - matriz de momento MLS;
  - multiplicadores de Lagrange;
  - GMRES;
  - distribuição regular versus não uniforme de nós.



## Importante - avaliado

[x] Para validação contra a solução analítica seno-seno-sinh do problema do cubo,
as arestas superiores devem herdar a condição das paredes laterais, V = 0.
A condição V = V0 é aplicada somente no interior aberto da face z = L.
Essa política evita impor, no mesmo ponto geométrico, duas condições de
Dirichlet incompatíveis e torna a fronteira numérica compatível com a
série analítica usada como referência.

Concluído em 2026-05-06. A política passou a ser global nos geradores
Dirichlet regular e por nuvem arbitrária. Com `15x15x15`, o erro global caiu
de `5.530370e-01` para `2.365161e-02` (erro relativo global), e o erro máximo
absoluto 3D caiu de `~V0` para `4.642650e-01`, mantendo `support_lt_4 = 0`,
`mls_failures = 0` e `gmres_converged = YES`.
