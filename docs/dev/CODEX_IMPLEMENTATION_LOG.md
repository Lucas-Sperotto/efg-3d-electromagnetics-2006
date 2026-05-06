# CODEX_IMPLEMENTATION_LOG — Registro de implementação

## Data

2026-05-05 22:19:00 -03

## Objetivo da rodada

Adicionar e executar os refinamentos regulares 13x13x13 e 15x15x15 no app
`reproduce_cube_sparse`, mantendo os casos 5x5x5 e 11x11x11, sem adicionar
nuvem não uniforme, precondicionador, octree, exportação de plano ou gráficos.

## Arquivos modificados

- `apps/reproduce_cube_sparse.c`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `docs/dev/TEST_REPORT.md`
- `TODO.md`

## Arquivos criados

- Nenhum arquivo versionado.

## Saída numérica CSV

O app grava o resumo numérico em:

```text
data/output/reproduce_cube_sparse_report.csv
```

O arquivo é ignorado por `.gitignore`, mas foi gerado na execução padrão.

## Resumo técnico das alterações

- Adicionados seletores de CLI:
  - `--case refine13`
  - `--case refine15`
- Atualizado `--case all` e a execução sem argumentos para rodar:
  - `sanity` 5x5x5;
  - `target` 11x11x11;
  - `refine13` 13x13x13;
  - `refine15` 15x15x15.
- Mantido `run_dense_cmp=1` apenas no caso 5x5x5.
- Mantido `run_dense_cmp=0` nos casos 11x11x11, 13x13x13 e 15x15x15.
- Configurados os refinamentos 13x13x13 e 15x15x15 com:
  - `restart=300`;
  - `max_iter=20000`;
  - `gmres_tol=1e-9`;
  - integração 15x15x15;
  - quadratura Gauss 2x2x2.
- Adicionado `max_cond_estimate` ao relatório obrigatório e ao CSV.

## Comandos executados

```bash
cmake --build build
/usr/bin/ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse --help
./build/reproduce_cube_sparse --case refine13
./build/reproduce_cube_sparse --case refine15
./build/reproduce_cube_sparse
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
Total Test time (real) = 0.14 sec
```

## Resultado refine13

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

## Resultado refine15

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

## Comparação 11/13/15

```text
11x11x11: interior error = 5.552359e-02, GMRES iters = 67, max_cond = 6.760583e+02
13x13x13: interior error = 4.566296e-02, GMRES iters = 66, max_cond = 9.004260e+02
15x15x15: interior error = 2.812108e-02, GMRES iters = 70, max_cond = 8.204838e+02
```

O erro interior não aumentou; caiu de ~5.55% para ~4.57% e depois para ~2.81%.

## Problemas encontrados

- Nenhum bloqueio numérico nos refinamentos:
  - `support_lt_4 = 0`;
  - `support_lt_8 = 0`;
  - `mls_failures = 0`;
  - `gmres_converged = YES`.
- O `ctest` sem caminho continua dependente do ambiente local; a validação foi
  feita com `/usr/bin/ctest`, conforme já documentado.

## Decisões tomadas

- Não foi implementada montagem direta esparsa de `K`; o app continua usando
  o caminho didático já validado.
- Não foi adicionado precondicionador porque GMRES convergiu sem ele.
- Não foi exportado o plano `x = 5.33`; isso fica para a próxima etapa.

## Próxima etapa recomendada

Exportar o plano `x = 5.33`, com CSV contendo `y,z,V_num,V_exact,abs_error`,
pois os casos regulares 11x11x11, 13x13x13 e 15x15x15 passaram nos critérios
de suporte, MLS e GMRES.
