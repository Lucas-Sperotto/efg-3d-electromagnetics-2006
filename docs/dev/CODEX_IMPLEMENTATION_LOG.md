# CODEX_IMPLEMENTATION_LOG — Registro de implementação

## Data

Preencher ao executar.

## Objetivo da rodada

Implementar ou ajustar o pipeline esparso para o problema regular do cubo eletrostático 3D.

## Arquivos modificados

- Preencher.

## Arquivos criados

- Preencher.

## Resumo técnico das alterações

Preencher com itens objetivos.

Exemplo:

- Criado app `reproduce_cube_sparse`.
- Adicionado alvo CMake `reproduce_cube_sparse`.
- Integrada montagem esparsa de `K`.
- Integrada montagem esparsa do sistema aumentado.
- Integrado solver GMRES.
- Adicionado relatório de erro e diagnóstico MLS.

## Comandos executados

```bash
cmake --build build
ctest --test-dir build --output-on-failure
./build/reproduce_cube_sparse
```

## Resultado do build

```text
Preencher.
```

## Resultado dos testes

```text
Preencher.
```

## Resultado do caso pequeno

Configuração:

```text
nodes:
cells:
constraints:
```

Saída:

```text
Preencher.
```

## Resultado do caso 11x11x11

```text
nodes:
constraints:
augmented_size:
K_nnz:
A_aug_nnz:
support_min:
support_mean:
support_max:
support_lt_4:
support_lt_8:
mls_failures:
gmres_converged:
gmres_iterations:
residual_initial:
residual_final:
relative_error_global:
relative_error_internal:
max_abs_error:
assembly_time_s:
solve_time_s:
```

## Problemas encontrados

* Preencher.

## Decisões tomadas

* Preencher.

## Pendências deixadas para a próxima IA

* Preencher.

## Recomendação para o Gemini

Revisar:

* clareza do app;
* mensagens de relatório;
* comentários didáticos;
* consistência do TODO;
* documentação dos resultados.
