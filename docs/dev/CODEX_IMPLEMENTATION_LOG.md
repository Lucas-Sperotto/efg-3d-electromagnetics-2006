# CODEX_IMPLEMENTATION_LOG — Registro de implementação

## Data

2026-05-05 21:55:38 -03

## Objetivo da rodada

Implementar o próximo passo da auditoria do Claude para o app esparso do cubo
eletrostático regular: manter o pipeline sparse/GMRES, executar os casos
5x5x5 e 11x11x11, imprimir o relatório obrigatório e registrar os resultados.

## Arquivos modificados

- `apps/reproduce_cube_sparse.c`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `docs/dev/TEST_REPORT.md`

## Arquivos criados

- Nenhum.

## Resumo técnico das alterações

- Mantido o alvo CMake existente `reproduce_cube_sparse`.
- Adicionada CLI mínima ao app: `--case sanity`, `--case target`, `--case all`,
  além de `--help`.
- Ajustada a estimativa inicial de capacidade de `K_coo` de `30 * node_count`
  para `100 * node_count + 16`, conforme recomendação da auditoria.
- Adicionado bloco terminal `--- Required report ---` com os campos pedidos:
  `nodes`, `constraints`, `augmented_size`, `K_nnz`, `A_aug_nnz`,
  `support_min`, `support_mean`, `support_max`, `support_lt_4`,
  `support_lt_8`, `mls_failures`, `gmres_converged`, `gmres_iterations`,
  `residual_initial`, `residual_final`, `relative_error_global`,
  `relative_error_internal`, `max_abs_error`, `assembly_time_s`,
  `solve_time_s`.
- Adicionado CSV de resumo em `data/output/reproduce_cube_sparse_report.csv`
  com os mesmos campos do relatório obrigatório.
- Adicionados critérios de parada explícitos no app:
  `support_lt_4 > 0`, `mls_failures > 0` ou `gmres_converged = NO`.
- Mantida a versão densa e a comparação densa apenas no caso pequeno.
- Não foram implementados nuvem não uniforme, precondicionador, octree,
  montagem de gráficos ou exportação de plano.

## Comandos executados

```bash
cmake --build build
ctest --test-dir build --output-on-failure
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --help
./build/reproduce_cube_sparse --case sanity
./build/reproduce_cube_sparse --case target
./build/reproduce_cube_sparse
```

## Resultado do build

```text
[100%] Built target test_gmres
```

Build completo concluído sem erros.

## Resultado dos testes

O comando esperado `ctest --test-dir build --output-on-failure` encontrou um
problema de ambiente:

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

Foi usado o binário real do CMake:

```bash
/usr/bin/ctest --test-dir build --output-on-failure
```

Resultado:

```text
100% tests passed, 0 tests failed out of 28
Total Test time (real) = 0.14 sec
```

## Resultado do caso pequeno

Configuração:

```text
nodes: 5x5x5
cells: 5x5x5
quadrature: Gauss 2x2x2 por célula
solver: GMRES, com comparação contra solver denso
```

Saída obrigatória:

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

## Resultado do caso 11x11x11

Configuração:

```text
L: 10
V0: 10
nodes: 11x11x11
cells: 15x15x15
quadrature: Gauss 2x2x2 por célula
support: terceiro vizinho multiplicado por 1.5
boundary: Dirichlet via multiplicadores de Lagrange
solver: GMRES sem precondicionador
```

Saída obrigatória:

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

## CSV gerado

Arquivo:

```text
data/output/reproduce_cube_sparse_report.csv
```

Conteúdo após `./build/reproduce_cube_sparse`:

```text
label,nodes,constraints,augmented_size,K_nnz,A_aug_nnz,support_min,support_mean,support_max,support_lt_4,support_lt_8,mls_failures,gmres_converged,gmres_iterations,residual_initial,residual_final,relative_error_global,relative_error_internal,max_abs_error,assembly_time_s,solve_time_s
"5x5x5 sanity check (dense comparison enabled)",125,98,223,6859,9715,8,17.575999999999997,27,0,0,0,YES,89,50,4.502731637452715e-08,0.65566020220249721,0.2178038587457618,10.002496316879469,0.0099749999999999995,0.000602
"11x11x11 target (15x15x15 integration cells, GMRES only)",1331,602,1933,117649,137383,8,23.557629629629595,27,0,0,0,YES,67,110,8.698755847818432e-08,0.55330567475497328,0.055523591643745318,10.000000001957561,0.41241299999999997,0.0083280000000000003
```

## Problemas encontrados

- O comando `ctest` sem caminho aponta para `/home/sperotto/.local/bin/ctest`,
  que falha por ausência do módulo Python `cmake`.
- O binário `/usr/bin/ctest` executa a suíte corretamente.
- Nenhum bloqueio numérico foi encontrado nos casos executados:
  `support_lt_4 = 0`, `mls_failures = 0`, `gmres_converged = YES`.

## Decisões tomadas

- O app continua usando `K_dense` como intermediário antes de montar `K_coo`,
  conforme aceito pela auditoria para o alvo 11x11x11.
- `support_lt_8` foi reportado como pontos com suporte ativo menor que 8,
  incluindo eventuais pontos inválidos. Nos casos executados o valor foi zero.
- O caso 11x11x11 manteve `restart=200` e `max_iter=10000`, pois é o alvo
  inicial pedido nesta rodada e convergiu em 67 iterações.

## Pendências deixadas para a próxima IA

- Se o objetivo passar para 13x13x13 ou 15x15x15, ajustar os parâmetros GMRES
  recomendados pela auditoria para `restart=300` e `max_iter=20000`.
- Avaliar montagem esparsa direta de `K` sem intermediário denso apenas depois
  dos próximos refinamentos exigirem isso.
- Exportar o plano `x = 5.33` somente após revisão dos resultados atuais.

## Recomendação para o Gemini

Revisar:

- consistência do bloco `--- Required report ---`;
- critério de parada para `support_lt_4`, `mls_failures` e GMRES;
- decisão de manter 11x11x11 como caso-alvo inicial;
- documentação do problema de ambiente com `ctest`;
- se o próximo passo deve ser 13x13x13/15x15x15 ou exportação do plano
  `x = 5.33`.
