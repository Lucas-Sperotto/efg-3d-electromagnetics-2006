# CODEX_IMPLEMENTATION_LOG — Registro de implementação

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
