# CLAUDE_AUDIT — Auditoria matemática e arquitetural

## Data

Preencher ao executar.

## Objetivo da auditoria

Avaliar se o projeto está pronto para executar o caso regular do cubo com:

```text
11x11x11 nós
15x15x15 células de integração
Dirichlet via Lagrange
matriz esparsa
GMRES
```

## Arquivos analisados

* [ ] TODO.md
* [ ] README.md
* [ ] CMakeLists.txt
* [ ] include/
* [ ] src/
* [ ] apps/
* [ ] tests/
* [ ] docs/

## Pontos verificados

### 1. Diagnóstico MLS

* [ ] Conta nós ativos no suporte.
* [ ] Verifica cardinalidade mínima da base linear 3D.
* [ ] Registra `active_nodes < 4`.
* [ ] Registra `active_nodes < 8`.
* [ ] Registra falhas MLS.
* [ ] Registra pivôs pequenos ou estimativa de condicionamento.

Observações:

```text
Preencher.
```

### 2. Base MLS linear 3D

Base esperada:

```text
p(x,y,z) = [1, x, y, z]
```

Cardinalidade mínima teórica:

```text
4 nós ativos
```

Observações:

```text
Preencher.
```

### 3. Sistema aumentado

Forma esperada:

```text
[ K  G^T ]
[ G   0  ]
```

Verificações:

* [ ] Dimensões coerentes.
* [ ] Sinais coerentes.
* [ ] Bloco `G` coerente com restrições de Dirichlet.
* [ ] Bloco `G^T` coerente.
* [ ] Vetor do lado direito coerente.

Observações:

```text
Preencher.
```

### 4. GMRES

* [ ] Usa matriz CSR.
* [ ] Retorna número de iterações.
* [ ] Retorna resíduo inicial.
* [ ] Retorna resíduo final.
* [ ] Indica convergência ou falha.
* [ ] Possui tolerância explícita.

Observações:

```text
Preencher.
```

### 5. Comparação com solução analítica

* [ ] Erro global.
* [ ] Erro interno.
* [ ] Erro máximo absoluto.
* [ ] Tratamento claro de pontos de contorno.
* [ ] Domínio físico `L = 10`.
* [ ] Potencial superior `V0 = 10`.

Observações:

```text
Preencher.
```

## Riscos técnicos encontrados

1. Preencher.
2. Preencher.
3. Preencher.

## Critérios de aceitação para o Codex

O Codex só deve considerar a etapa concluída se:

```text
support_lt_4 = 0
mls_failures = 0
gmres_converged = true
residual_final aceitável
erro interno reportado
testes passando
```

## Tarefas recomendadas ao Codex

1. Criar ou ajustar `apps/reproduce_cube_sparse.c`.
2. Adicionar alvo `reproduce_cube_sparse` ao CMake.
3. Rodar caso pequeno de sanidade.
4. Rodar caso `11x11x11`.
5. Registrar resultados em `docs/dev/TEST_REPORT.md`.
6. Registrar alterações em `docs/dev/CODEX_IMPLEMENTATION_LOG.md`.

## O que não deve ser feito ainda

* Nuvens não uniformes.
* Octree.
* Precondicionador.
* Gráficos finais.
* Refatoração ampla.
* Mudança de APIs públicas sem necessidade.

## Recomendação

```text
AVANÇAR / AVANÇAR COM CUIDADO / NÃO AVANÇAR
```
